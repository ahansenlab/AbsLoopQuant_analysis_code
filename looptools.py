import os
import numpy as np
import pandas as pd
import cv2
from scipy.ndimage import grey_dilation
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Circle
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rcParams
rcParams['font.family'] = 'Arial'
rcParams['font.size'] = 8
import cooltools
import warnings
from cooltools.lib import plotting

def get_image(clr, chr_name, left, right, pad, balance=True):
    """Get the image of the loop +/- the provided pad distance."""
    return clr.matrix(balance=balance).fetch(f'{chr_name}:{left-pad}-{left+pad}',f'{chr_name}:{right-pad}-{right+pad}')

def get_s_px_matrix(loop_size, pad, res):
    """Get a matrix of s (genomic separation) in units of the cooler resolution, whose support corresponds to the image."""
    loop_size_px = loop_size // res
    a = pad // res
    y_px, x_px = np.meshgrid(np.arange(-a, a+1),np.arange(-a, a+1))
    s_px_matrix = loop_size_px+y_px-x_px  # genomic separation in units of res
    s_px_matrix[s_px_matrix<0] = 0  # don't allow negative values of s (set them to 0)
    return s_px_matrix

def bin_scatter_plot(x, y, nbins=10):
    """Given x and y data, split the x-axis into nbins bins, and return points whose x and y values are the averages of the x and y values of the data points in the bins. Also return the y-value's standard error (stdev/sqrt(N))."""
    epsilon = (np.max(x)-np.min(x))*0.01
    bin_edges = np.linspace(np.min(x), np.max(x)+epsilon, nbins+1)[1:]
    bin_num_x = np.digitize(x, bin_edges)
    x_binned = np.zeros(nbins)
    y_binned = np.zeros(nbins)
    y_err_binned = np.zeros(nbins)
    for bin_num in range(0, nbins):
        x_binned[bin_num] = np.mean(x[bin_num_x==bin_num])
        y_binned[bin_num] = np.mean(y[bin_num_x==bin_num])
        y_err_binned[bin_num] = np.std(y[bin_num_x==bin_num], ddof=1)/np.sqrt(np.sum(bin_num_x==bin_num))
    return x_binned, y_binned, y_err_binned

def calculate_and_save_avg_Ps_curve(clr, nproc=1, max_sep=3000000, output_filename=None):
    """Calculate and save the P(s) curve (average across chromosomes, weighted by chromosome size) and save to file. max_sep is the maximum genomic separation in bp to which to calculate the P(s) curve."""

    res = clr.binsize
    
    # calculate the P(s) curve
    cvd = cooltools.expected_cis(clr=clr, smooth=True, aggregate_smoothed=True, nproc=nproc)
    cvd['s_bp'] = cvd['dist'] * res
    
    chr_names = [chrom_name for chrom_name in clr.chromnames if len(chrom_name)>3 and (chrom_name[3:].isnumeric() or chrom_name[3:]=='X')]  # only take numbered choromosomes and chrX

    # average across chromosomes
    P_s_data_all_chrs = np.zeros((1+max_sep//res, len(chr_names)))
    for i, chr_name in enumerate(chr_names):
        P_s_data_chr = cvd.loc[cvd['region1']==chr_name,np.array(['s_bp','balanced.avg'])]
        P_s_data_chr = P_s_data_chr.loc[P_s_data_chr['s_bp']<=max_sep]
        P_s_data_chr = P_s_data_chr.sort_values('s_bp')
        assert(np.all(P_s_data_chr['s_bp']==np.arange(0,max_sep+1,res)))
        P_s_data_all_chrs[:,i] = P_s_data_chr['balanced.avg']

    chrom_weights = clr.chromsizes[chr_names].values
    chrom_weights = chrom_weights/np.sum(chrom_weights)
    
    P_s_data_averaged_chrs = np.average(P_s_data_all_chrs, axis=1, weights=chrom_weights)

    # save as txt file
    if output_filename is None:
        output_filename = f"P_s_{clr.filename.split('/')[-1].split('.')[0]}_{res}bp.txt"
    np.savetxt(output_filename, P_s_data_averaged_chrs)

class LoopQuantifier:
    def __init__(self, clr, P_s_values, gaussian_blur_sigma_px=10, outlier_removal_radius_px=10, ignore_diag_cutoff_px=5, na_stripe_dist_to_center_px_cutoff=5, footprint=np.array([[0,1,1,1,0],[1,1,1,1,1],[1,1,1,1,1],[1,1,1,1,1],[0,1,1,1,0]])):
        self.clr = clr
        self.res = clr.binsize  # resolution in bp
        self.P_s_values = P_s_values
        self.Ps_pd = pd.Series(P_s_values, index=np.arange(len(P_s_values))*self.res)  # make pandas version to allow indexing with real-valued numbers (in curve_fit)
        self.gaussian_blur_sigma_px = gaussian_blur_sigma_px
        self.outlier_removal_radius_px = outlier_removal_radius_px
        self.ignore_diag_cutoff_px = ignore_diag_cutoff_px
        self.na_stripe_dist_to_center_px_cutoff = na_stripe_dist_to_center_px_cutoff
        self.footprint = footprint  # for finding local maxima
        
        self.local_region_size = None  # initialize as None, will be assigned during loop quantification
        self.quant_region_size = None  # initialize as None, will be assigned during loop quantification
    
    def snap_coordinate_to_bin_center(self, x):
        """Given a genomic coordinate, "snap" it to the center of the bin it is located in."""
        return x//self.res*self.res + self.res//2
        
    def generate_precomputed_matrices(self, local_region_size, quant_region_size):
        """
        Pre-compute matrices necessary for loop quantification that depend on local_region_size and quant_region_size, most importantly:
        (1) x_px and y_px: matrices of Cartesian coordinates of pixels centered at (0, 0) that is the same size as the local region
        (2) boolean matrix describing the diamond mask for the local region
        (3) boolean matrix describing the circular mask for the quantification region
        and define constant values necessary for loop quantification.
        """
        a = local_region_size//self.res
        self.y_px, self.x_px = np.meshgrid(np.arange(-a, a+1),np.arange(-a, a+1))
        self.diamond = np.abs(self.x_px)+np.abs(self.y_px)<=a
        self.diamond_vertices = [(0,a),(a,0),(0,-a),(-a,0)]
        self.circle = np.sqrt(self.x_px**2+self.y_px**2)<=quant_region_size//self.res
        self.circle_expanded = np.sqrt(self.x_px**2+self.y_px**2)<=quant_region_size//self.res+1  # slightly larger, for plotting purposes
        self.extent_px = np.array((-a-0.5, a+0.5, a+0.5, -a-0.5))
        
        self.local_region_size = local_region_size
        self.quant_region_size = quant_region_size
    
    def resolve_NAs(self, mat):
        """Resolve NAs in a matrix by replacing them by the median value."""
        matrix_dimension = mat.shape[0]  # height/width of square image
        ver_na_stripe_indices = np.where(np.sum(np.isnan(mat),0)==matrix_dimension)[0]  # get indices of NA stripes
        hor_na_stripe_indices = np.where(np.sum(np.isnan(mat),1)==matrix_dimension)[0]
        any_na_stripes = len(ver_na_stripe_indices)>0 or len(hor_na_stripe_indices)>0  # boolean indicating whether or not there are any NA stripes

        if any_na_stripes:
            middle_index = matrix_dimension//2
            ver_na_stripe_indices_from_middle = np.abs(ver_na_stripe_indices - middle_index)
            hor_na_stripe_indices_from_middle = np.abs(hor_na_stripe_indices - middle_index)
            if np.any(ver_na_stripe_indices_from_middle<=self.na_stripe_dist_to_center_px_cutoff) or np.any(hor_na_stripe_indices_from_middle<=self.na_stripe_dist_to_center_px_cutoff):
                warnings.warn("Removal of NaN values failed; NaN-valued pixels too close to center!", stacklevel=2)
                return mat*np.nan  # matrix of all NA's
            mat_NAs_removed = np.nan_to_num(mat, nan=np.nanmedian(mat[self.diamond]))  # replace NA values with median within the diamond
        else:
            mat_NAs_removed = mat
            
        return mat_NAs_removed
    
    def detect_outliers(self, chr_name, left, right, local_region_size=50000, quant_region_size=10000, clr_for_outlier_detection=None, P_s_values_for_outlier_detection=None, k_min=None):
        
        # default values
        if clr_for_outlier_detection is None:
            clr_for_outlier_detection = self.clr
        if P_s_values_for_outlier_detection is None:
            P_s_values_for_outlier_detection = self.P_s_values
        if k_min is None:
            k_min = 2
        
        # if the stored values of local_region_size and quant_region_size do not match the requested values,
        # regenerate the precomputed matrices and store the new values
        if (self.local_region_size!=local_region_size) or (self.quant_region_size!=quant_region_size):
            self.generate_precomputed_matrices(local_region_size, quant_region_size)
        
        left = self.snap_coordinate_to_bin_center(left)
        right = self.snap_coordinate_to_bin_center(right)
        
        # get image and expected background image
        img = get_image(clr_for_outlier_detection, chr_name, left, right, pad=self.local_region_size)
        s_px_matrix = get_s_px_matrix(loop_size=right-left, pad=self.local_region_size, res=self.res)
        bg_img = P_s_values_for_outlier_detection[s_px_matrix]
        
        # make all pixels near diagonal NA
        img[s_px_matrix<=self.ignore_diag_cutoff_px] = np.nan
        bg_img[s_px_matrix<=self.ignore_diag_cutoff_px] = np.nan
        
        # divide the image by the expected background
        img_over_bg = img/bg_img
        
        # resolve NAs
        img_over_bg_NAs_removed = self.resolve_NAs(img_over_bg)
            
        # blur image
        ksize = int(np.ceil(3*self.gaussian_blur_sigma_px)//2*2+1)  # round up to next odd integer >= 3 sigma
        img_over_bg_blurred = cv2.GaussianBlur(img_over_bg_NAs_removed,ksize=(ksize,ksize),sigmaX=self.gaussian_blur_sigma_px)
        
        # crop blurred image to diamond before finding local maxima
        img_over_bg_blurred[~self.diamond] = np.nan
        
        # find local maxima
        local_maxima_bool = img_over_bg_blurred==grey_dilation(img_over_bg_blurred, footprint=self.footprint)
        strong_local_maxima_bool = np.logical_and(local_maxima_bool, img_over_bg_blurred>k_min*np.nanmedian(img_over_bg_blurred))
        
        return strong_local_maxima_bool
    
    def quantify_loop(self, chr_name, left, right, coords_convert_function=None, local_region_size=50000, quant_region_size=10000, clr_for_outlier_detection=None, P_s_values_for_outlier_detection=None, convert_coords_outliers=False, k_min=None, outliers_to_remove=None, return_all=False, show_plot=False):
        
        # if the stored values of local_region_size and quant_region_size do not match the requested values,
        # regenerate the precomputed matrices and store the new values
        if (self.local_region_size!=local_region_size) or (self.quant_region_size!=quant_region_size):
            self.generate_precomputed_matrices(local_region_size, quant_region_size)
        
        if coords_convert_function is None:
            left_quant = left
            right_quant = right
        else:
            left_quant = coords_convert_function(chr_name, left)
            right_quant = coords_convert_function(chr_name, right)
            
        if convert_coords_outliers:
            left_outliers = left_quant
            right_outliers = right_quant
        else:
            left_outliers = left
            right_outliers = right
        
        left_quant = self.snap_coordinate_to_bin_center(left_quant)
        right_quant = self.snap_coordinate_to_bin_center(right_quant)
        left_outliers = self.snap_coordinate_to_bin_center(left_outliers)
        right_outliers = self.snap_coordinate_to_bin_center(right_outliers)
        
        # get image and expected background image
        img = get_image(self.clr, chr_name, left_quant, right_quant, pad=self.local_region_size)
        s_px_matrix = get_s_px_matrix(loop_size=right_quant-left_quant, pad=self.local_region_size, res=self.res)
        bg_img = self.P_s_values[s_px_matrix]
        
         # divide the image by the expected background
        img_over_bg = img/bg_img
        
        # resolve NAs
        img_NAs_removed = self.resolve_NAs(img_over_bg) * bg_img
        
        # crop image to diamond
        img[~self.diamond] = np.nan
        img_NAs_removed[~self.diamond] = np.nan
        
        # remove outliers
        img_outliers_removed = img_NAs_removed.copy()
        if outliers_to_remove is None:
            outliers_to_remove = np.where(self.detect_outliers(chr_name, left_outliers, right_outliers, local_region_size, quant_region_size, clr_for_outlier_detection, P_s_values_for_outlier_detection, k_min))
        i_local_max_arr, j_local_max_arr = outliers_to_remove
        for k in range(len(i_local_max_arr)):
            x_local_max = self.x_px[:,0][i_local_max_arr[k]]
            y_local_max = self.y_px[0,:][j_local_max_arr[k]]
            dist_to_local_max = np.sqrt((self.x_px-x_local_max)**2+(self.y_px-y_local_max)**2)
            img_outliers_removed[dist_to_local_max<=self.outlier_removal_radius_px] = np.nan
            
        # fit P(s) curve
        s_to_fit = s_px_matrix[s_px_matrix>self.ignore_diag_cutoff_px].flatten() * self.res
        P_s_to_fit = img_outliers_removed[s_px_matrix>self.ignore_diag_cutoff_px].flatten()
        s_to_fit, P_s_to_fit = s_to_fit[np.logical_and(~np.isnan(s_to_fit), ~np.isnan(P_s_to_fit))], P_s_to_fit[np.logical_and(~np.isnan(s_to_fit), ~np.isnan(P_s_to_fit))]
        c_best_fit = curve_fit(lambda s,c: c*self.Ps_pd[s], s_to_fit, P_s_to_fit)[0][0]
        
        # get local background
        local_bg_img = bg_img * c_best_fit
             
        # subtract local background from image
        img_local_bg_subtracted = img_NAs_removed - local_bg_img
        
        # crop to circle
        img_local_bg_subtracted[~self.circle_expanded] = np.nan
        loop_quantification_score = np.sum(img_local_bg_subtracted[self.circle])
              
        if show_plot:
            self.plot_quantification(chr_name, left, right, img, img_over_bg, img_outliers_removed, img_local_bg_subtracted, s_to_fit, P_s_to_fit, c_best_fit, coords_convert_function)
            
        if return_all:
            output = {'loop_quantification_score':loop_quantification_score,
                      'img':img,
                      'img_over_bg':img_over_bg,
                      'img_NAs_removed':img_NAs_removed,
                      'img_outliers_removed':img_outliers_removed,
                      'local_bg_img':local_bg_img,
                      'img_local_bg_subtracted':img_local_bg_subtracted,
                      's_to_fit':s_to_fit,
                      'P_s_to_fit':P_s_to_fit,
                      'c_best_fit':c_best_fit}
            return output
            
        return loop_quantification_score
    
    def plot_quantification(self, chr_name, left, right, img, img_over_bg, img_outliers_removed, img_local_bg_subtracted, s_to_fit, P_s_to_fit, c_best_fit, coords_convert_function=None):
        
        if coords_convert_function is not None:  # for first plot only
            left_large_plot = coords_convert_function(chr_name, left)
            right_large_plot = coords_convert_function(chr_name, right)
        else:
            left_large_plot = left
            right_large_plot = right
        
        # make figure and subplots
        fig, axs = plt.subplots(3, 2, figsize=(6, 8.5))
        ax1, ax2, ax3, ax4, ax5, ax6 = axs.flatten()
       
        # position colorbar axes
        pos2 = ax2.get_position()
        cax2 = fig.add_axes([pos2.xmax-(pos2.xmax-pos2.xmin)*0.12, pos2.ymin, (pos2.xmax-pos2.xmin)*0.06, (pos2.ymax-pos2.ymin)*0.3])        
        pos3 = ax3.get_position()
        cax3 = fig.add_axes([pos3.xmax-(pos3.xmax-pos3.xmin)*0.12, pos3.ymin, (pos3.xmax-pos3.xmin)*0.06, (pos3.ymax-pos3.ymin)*0.3])
        pos4 = ax4.get_position()
        cax4 = fig.add_axes([pos4.xmax-(pos4.xmax-pos4.xmin)*0.12, pos4.ymin, (pos4.xmax-pos4.xmin)*0.06, (pos4.ymax-pos4.ymin)*0.3])
        pos6 = ax6.get_position()
        cax6 = fig.add_axes([pos6.xmax-(pos6.xmax-pos6.xmin)*0.12, pos6.ymin, (pos6.xmax-pos6.xmin)*0.06, (pos6.ymax-pos6.ymin)*0.3])
      
        # plot images
        im2 = ax2.imshow(img, extent=self.extent_px, cmap='fall', vmin=0, vmax=np.nanquantile(img, 0.99))
        fig.colorbar(im2, cax=cax2, orientation='vertical')
        ax2.set_title(f'local region (loop $\pm$ {int(self.local_region_size/1e3)} kb)')
        
        eps = int((right_large_plot-left_large_plot)*0.1)
        zoomed_cmap = self.clr.matrix().fetch(f'{chr_name}:{left_large_plot-eps}-{right_large_plot+eps}',f'{chr_name}:{left_large_plot-eps}-{right_large_plot+eps}')
        im1 = ax1.imshow(zoomed_cmap, extent=np.array([right_large_plot+eps,left_large_plot-eps,left_large_plot-eps,right_large_plot+eps])/1e6, origin='upper', cmap='fall', vmin=0, vmax=np.nanquantile(img, 0.99))
        ax1.set_title('loop and anchors')
        
        im3 = ax3.imshow(img_over_bg, extent=self.extent_px, cmap='fall', vmin=0, vmax=np.nanquantile(img_over_bg, 0.99))
        fig.colorbar(im3, cax=cax3, orientation='vertical')
        ax3.set_title('image divided by background')
        
        im4 = ax4.imshow(img_outliers_removed, extent=self.extent_px, cmap='fall', vmin=0, vmax=np.nanquantile(img, 0.99))
        fig.colorbar(im4, cax=cax4, orientation='vertical')
        ax4.set_title('local region, outliers removed')
        
        for ax, im in zip([ax2, ax3, ax4], [im2, im3, im4]):
            ax.add_patch(Polygon(self.diamond_vertices, edgecolor='black', linewidth=1, fill=None))
            im.set_clip_path(Polygon(self.diamond_vertices, transform=ax.transData))
            ax.axis('off')
            
        x_binned, y_binned, y_std_binned = bin_scatter_plot(s_to_fit, P_s_to_fit, nbins=20)
        ax5.errorbar(x_binned/1e6, y_binned, y_std_binned, color='#292929', fmt='o', ms=3, linewidth=1, capsize=2)
        x_best_fit = np.arange(np.min(s_to_fit), np.max(s_to_fit)+1, self.res)
        y_best_fit = self.Ps_pd[x_best_fit] * c_best_fit
        ax5.plot(x_best_fit/1e6, y_best_fit, color='#2222ee')
        ax5.set_title('P(s) curve fitting')
            
        im6 = ax6.imshow(img_local_bg_subtracted, extent=self.extent_px, cmap='bwr', vmax=np.nanmax(np.abs(img_local_bg_subtracted)), vmin=-np.nanmax(np.abs(img_local_bg_subtracted)))
        ax6.add_patch(Circle((0, 0), radius=self.quant_region_size/self.res, edgecolor='black', linewidth=1, fill=None))
        im6.set_clip_path(Circle((0, 0), radius=self.quant_region_size/self.res, transform=ax6.transData))
        ax6.set_xlim(-self.quant_region_size/self.res-1,self.quant_region_size/self.res+1)
        ax6.set_ylim(self.quant_region_size/self.res+1, -self.quant_region_size/self.res-1)
        fig.colorbar(im6, cax=cax6, orientation='vertical')
        ax6.axis('off')
        ax6.set_title(f'quantification region ({int(self.quant_region_size/1e3)} kb radius)')
        
        plt.suptitle(f'{chr_name}:{left/1e6:.3f}-{right/1e6:.3f} Mb', fontsize=12)
        
        plt.subplots_adjust(hspace=0.3)