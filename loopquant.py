class LoopQuantifier:
    def __init__(self, clr, P_s_values, gaussian_blur_sigma_px=10, outlier_removal_radius_px=10, ignore_diag_cutoff_px=5, na_stripe_dist_to_center_px_cutoff=5):
        self.clr = clr
        self.res = clr.binsize  # resolution in bp
        self.P_s_values = P_s_values
        self.Ps_pd = pd.Series(P_s_values, index=np.arange(len(P_s_values))*res)  # make pandas version to allow indexing with real-valued numbers (in curve_fit)
        self.gaussian_blur_sigma_px = gaussian_blur_sigma_px
        self.outlier_removal_radius_px = outlier_removal_radius_px
        self.ignore_diag_cutoff_px = ignore_diag_cutoff_px
        self.na_stripe_dist_to_center_px_cutoff = na_stripe_dist_to_center_px_cutoff
        self.local_region_size = None
        self.quant_region_size = None
        
        self.FOOTPRINT = np.array([[0,1,1,1,0],[1,1,1,1,1],[1,1,1,1,1],[1,1,1,1,1],[0,1,1,1,0]])  # for finding local maxima
    
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
    
    def get_image_local_region(self, chr_name, left, right, clr=self.clr):
        """Get the image of the local region surrounding the loop."""
        img = clr.matrix(balance=True).fetch(f'{chr_name}:{left-self.local_region_size}-{left+self.local_region_size}',f'{chr_name}:{right-self.local_region_size}-{right+self.local_region_size}')
        return img
    
    def get_s_px_matrix(self, chr_name, loop_size):
        """Get a matrix of s (genomic separation) in units of the cooler resolution, whose support corresponds to the local region."""
        loop_size_px = loop_size // self.res
        s_px_matrix = loop_size_px+self.y_px-self.x_px  # genomic separation in units of res
        s_px_matrix[s_px_matrix<0] = 0  # don't allow negative values of s (set them to 0)
        return s_px_matrix
    
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
                return [[],[]]
            mat_NAs_removed = np.nan_to_num(mat, nan=np.nanmedian(mat[self.diamond]))  # replace NA values with median within the diamond
        else:
            mat_NAs_removed = mat
            
        return mat_NAs_removed
        
    
    def detect_outliers(self, chr_name, left, right, local_region_size, quant_region_size, clr_for_outlier_detection=self.clr, P_s_values_for_outlier_detection=self.P_s_values):
        
        # if the stored values of local_region_size and quant_region_size do not match the requested values,
        # regenerate the precomputed matrices and store the new values
        if (self.local_region_size!=local_region_size) or (self.quant_region_size!=quant_region_size):
            self.generate_precomputed_matrices(local_region_size, quant_region_size)
        
        left = self.snap_coordinate_to_bin_center(left)
        right = self.snap_coordinate_to_bin_center(right)
        
        # get image and expected background image
        img = get_image_local_region(chr_name, left, right, clr=clr_for_outlier_detection)
        s_px_matrix = get_s_px_matrix(chr_name, loop_size=right-left)
        bg_img = P_s_values_for_outlier_detection[s_px_matrix]
        
        # make all pixels near diagonal NA
        img[s_px_matrix<=self.ignore_diag_cutoff_px] = np.nan
        bg_img[s_px_matrix<=self.ignore_diag_cutoff_px] = np.nan
        
        # divide the image by the expected global background
        img_over_bg = img/bg_img
        
        # resolve NAs
        img_over_bg_NAs_removed = resolve_NAs(img_over_bg)
            
        # blur image
        ksize = int(np.ceil(3*self.gaussian_blur_sigma_px)//2*2+1)  # round up to next odd integer >= 3 sigma
        img_over_bg_blurred = cv2.GaussianBlur(img_over_bg_NAs_removed,ksize=(ksize,ksize),sigmaX=self.gaussian_blur_sigma_px)
        
        # crop blurred image to diamond before finding local maxima
        img_over_bg_blurred[~self.diamond] = np.nan
        
        # find local maxima
        local_maxima_bool = img_over_bg_blurred==grey_dilation(img_over_bg_blurred, footprint=self.FOOTPRINT)
        strong_local_maxima_bool = np.logical_and(local_maxima_bool, img_over_bg_blurred>2*np.nanmedian(img_over_bg_blurred))
        
        # store variables for analysis
        self.last_clr_name_outlier_detection = os.path.basename(clr.filename).split('.')[0]
        self.last_chr_name = chr_name
        self.last_left = left
        self.last_right = right
        self.last_img = img
        self.last_bg_img = bg_img
        self.last_img_over_bg = img_over_bg
        self.last_img_over_bg_NAs_removed = img_over_bg_NAs_removed
        self.last_local_maxima_bool = local_maxima_bool
        self.last_img_over_bg_blurred = img_over_bg_blurred
        self.last_strong_local_maxima_bool = strong_local_maxima_bool
        
        return strong_local_maxima_bool
    
    def quantify_loop(self, clr, chr_name, left, right, outliers_to_remove):
        # adjust left and right to be in bin centers
        left = left//self.res*self.res + self.res//2
        right = right//self.res*self.res + self.res//2
        
        # get the image
        img = clr.matrix(balance=True).fetch(f'{chr_name}:{left-self.local_region_size}-{left+self.local_region_size}',f'{chr_name}:{right-self.local_region_size}-{right+self.local_region_size}')
        
        # get s_px matrix
        loop_size_px = right//self.res-left//self.res
        s_px_matrix = loop_size_px+self.y_px-self.x_px  # genomic separation in units of res
        s_px_matrix[s_px_matrix<0] = 0  # don't allow negative values of s
        
        # get the expected global background image
        bg_img = self.P_s_data[s_px_matrix]
        
        # in the image and background image, make all pixels near diagonal NA
        img[s_px_matrix<=self.ignore_diag_cutoff_px] = np.nan
        bg_img[s_px_matrix<=self.ignore_diag_cutoff_px] = np.nan
        
        # resolve NAs
        length_of_img = img.shape[0]  # height/width of square image
        ver_na_stripe_indices = np.where(np.sum(np.isnan(img),0)==length_of_img)[0]  # get indices of NA stripes
        hor_na_stripe_indices = np.where(np.sum(np.isnan(img),1)==length_of_img)[0]
        any_na_stripes = len(ver_na_stripe_indices)>0 or len(hor_na_stripe_indices)>0

        if any_na_stripes:
            middle_index = length_of_img//2
            ver_na_stripe_indices_from_middle = np.abs(ver_na_stripe_indices - middle_index)
            hor_na_stripe_indices_from_middle = np.abs(hor_na_stripe_indices - middle_index)
            if np.any(ver_na_stripe_indices_from_middle<=self.na_stripe_dist_to_center_px_cutoff) or np.any(hor_na_stripe_indices_from_middle<=self.na_stripe_dist_to_center_px_cutoff):
                warnings.warn("Outlier detection failed; NaN-valued pixels too close to center!", stacklevel=2)
                return None
            img_over_bg = img/bg_img
            img_over_bg_NAs_removed = np.nan_to_num(img_over_bg, nan=np.nanmedian(img_over_bg[self.diamond]))  # replace NA values with median within the diamond
            img_NAs_removed = img_over_bg_NAs_removed*bg_img
        else:
            img_NAs_removed = img
        
        # crop image to diamond
        img[~self.diamond] = np.nan
        img_NAs_removed[~self.diamond] = np.nan
        
        # remove outliers
        img_outliers_removed = img_NAs_removed.copy()
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
        
        # store variables for analysis
        self.last_clr_name_quantification = os.path.basename(clr.filename).split('.')[0]
        self.last_chr_name = chr_name
        self.last_left = left
        self.last_right = right
        self.last_img = img
        self.last_bg_img = bg_img
        self.last_img_NAs_removed = img_NAs_removed
        self.last_img_outliers_removed = img_outliers_removed
        self.last_s_to_fit = s_to_fit
        self.last_P_s_to_fit = P_s_to_fit
        self.last_c_best_fit = c_best_fit
        self.last_local_bg_img = local_bg_img
        self.last_img_local_bg_subtracted = img_local_bg_subtracted
        self.last_loop_quantification_score = loop_quantification_score
        
        return loop_quantification_score
    
    def plot_outlier_detection(self):
        # plot of raw image, img over background, blurred image, and blurred image with outliers
        # must be run immediately after outlier detection
        
        fig = plt.figure()
        gs = fig.add_gridspec(1, 4, hspace=0)
        ax1, ax2, ax3, ax4 = gs.subplots()
        im1 = ax1.imshow(self.last_img)
        im2 = ax2.imshow(self.last_img_over_bg)
        im3 = ax3.imshow(self.last_img_over_bg_blurred)
        im4 = ax4.imshow(self.last_img_over_bg_blurred)
        outlier_indices = np.where(self.last_strong_local_maxima_bool)
        ax4.scatter(outlier_indices[1], outlier_indices[0], facecolor='none', edgecolor='red', s=100)
    
    def plot_quantification(self):
        # plot of raw image, image with outliers removed, img with background subtracted, loop quantification window
        # must be run immediately after running quantification
        
        fig, axs = plt.subplots(2, 3, figsize=(8,5))
        axs = axs.flatten()
        ims = [None]*4

        ims[0] = axs[0].imshow(self.last_img, extent=self.extent_px, cmap='viridis', vmin=0, vmax=np.nanquantile(self.last_img, 1))
        ims[3] = axs[3].imshow(self.last_img_outliers_removed, extent=self.extent_px, cmap='viridis', vmin=0, vmax=np.nanquantile(self.last_img, 1))
        for i in [0,3]:
            axs[i].add_patch(Polygon(self.diamond_vertices, edgecolor='black', linewidth=1, fill=None))
            ims[i].set_clip_path(Polygon(self.diamond_vertices, transform=axs[i].transData))
            axs[i].axis('off')
        axs[0].set_title('local region (raw)')
        axs[3].set_title('local region (outliers removed)')

        x = self.last_s_to_fit
        y = self.last_P_s_to_fit
        axs[4].scatter(x, y, s=0.5, alpha=0.2, color='#e89c20')
        x_binned, y_binned, y_std_binned = bin_scatter_plot(x, y, nbins=20)
        axs[4].errorbar(x_binned, y_binned, y_std_binned, color='#292929', fmt='o', ms=3, linewidth=1, capsize=2)
        x_best_fit = np.arange(np.min(x), np.max(x)+1, self.res)
        y_best_fit = self.Ps_pd[x_best_fit] * self.last_c_best_fit
        axs[4].plot(x_best_fit, y_best_fit, color='#2222ee')
        axs[4].set_xlabel('s [bp]')
        axs[4].set_ylabel('P(s) [arb.]')
        axs[4].set_title('estimation of local background')
        axs[4].text(axs[4].get_xlim()[1], axs[4].get_ylim()[1]*0.95, f'c$=${self.last_c_best_fit:.2f}   ', ha='right', va='top')

        im_circle = axs[5].imshow(self.last_img_local_bg_subtracted, extent=self.extent_px, cmap='bwr', vmax=np.nanquantile(self.last_img_local_bg_subtracted, 1), vmin=-np.nanquantile(self.last_img_local_bg_subtracted, 1))
        axs[5].add_patch(Circle((0, 0), radius=self.quant_region_size/self.res, edgecolor='black', linewidth=1, fill=None))
        im_circle.set_clip_path(Circle((0, 0), radius=self.quant_region_size/self.res, transform=axs[5].transData))
        axs[5].set_xlim(-self.quant_region_size/self.res-1,self.quant_region_size/self.res+1)
        axs[5].set_ylim(self.quant_region_size/self.res+1, -self.quant_region_size/self.res-1)
        axs[5].axis('off')
        axs[5].set_title('loop quantification window\n(raw$-$local background)', pad=0)
        axs[5].text(0, self.quant_region_size/self.res+3, f'$\Sigma$ pixels$=${self.last_loop_quantification_score:.3f}', va='center', ha='center')

        plt.suptitle(f'{self.last_chr_name}:$\,${(self.last_left-self.res/2)/1e6:.3f}\u2013{(self.last_right-self.res/2)/1e6:.3f} Mb ({(self.last_right-self.last_left)/1e3:.0f} kb)\nsample: __, bin size: {self.res} bp')

        plt.tight_layout()