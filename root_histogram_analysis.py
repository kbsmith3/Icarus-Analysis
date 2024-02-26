import ROOT
import uproot
import numpy as np
import matplotlib.pyplot as plt
import os
import awkward 
import cython_function
import array

class Histogram_Python_Loader:
    def __init__(self, file_data_path = 'histogramsDATA.root', file_mc_path = 'histogramsMC.root', small_constant=1e-6):
        self.file_data_path = file_data_path
        self.file_mc_path = file_mc_path
        self.small_constant = small_constant
        self.ratio_histograms = {}
        self.histogram_data_arrays = {}
        self.histogram_mc_arrays = {}
        
        file = uproot.open(self.file_data_path)

        for key in file.keys():
            obj = file[key]

            if 'TH2' in obj.classname:  # Check for 2D histogram
                # Retrieve the histogram data and axis information
                hist_array, x_edges, y_edges = obj.to_numpy()
                x_title = obj.title.split(' ')[0]
                y_title = obj.title.split(' ')[1]
                hist_plane= key.split(' ')[1]
                hits = key.split(' ')[6]
                
                self.histogram_data_arrays[key] = {"title": obj.title, "histogram": hist_array, "x_edges": x_edges, "y_edges": y_edges}

        file = uproot.open(self.file_mc_path)

        for key in file.keys():
            obj = file[key]

            if 'TH2' in obj.classname:  # Check for 2D histogram
                # Retrieve the histogram data and axis information
                hist_array, x_edges, y_edges = obj.to_numpy()
                x_title = obj.title.split(' ')[0]
                y_title = obj.title.split(' ')[1]
                hist_plane= key.split(' ')[1]
                hits = key.split(' ')[6]
                
                self.histogram_mc_arrays[key] = {"title": obj.title, "histogram": hist_array, "x_edges": x_edges, "y_edges": y_edges}
                
class Combined_Histograms(Histogram_Python_Loader):
    def __init__(self, selected_x_title , selected_y_title ,selected_plane = False, min_data = 100):
        Histogram_Python_Loader.__init__(self)
        
        self.X = selected_x_title
        self.Y = selected_y_title
        
        self.min_data = min_data
        
        xedges = np.zeros(shape = (100))
        yedges = np.zeros(shape = (100))
        combined_array1 = np.zeros(shape = (100,100))
        combined_array2 = np.zeros(shape = (100,100))


        for key in self.histogram_data_arrays:

            x_title = self.histogram_data_arrays[key]["title"].split(' ')[0]
            y_title = self.histogram_data_arrays[key]["title"].split(' ')[1]
            hist_plane= key.split(' ')[1]
            hits = key.split(' ')[6]
            if selected_plane == False:
                plane = hist_plane
            else:
                plane = selected_plane

            # Check if the axes titles match the selected ones
            if x_title == selected_x_title and y_title == selected_y_title and (hist_plane == plane) and hits == "ALL":

                # Add the current histogram to the combined array for this pair of axes
                combined_array1 += self.histogram_data_arrays[key]["histogram"]
                xedges = self.histogram_data_arrays[key]['x_edges']
                yedges = self.histogram_data_arrays[key]['y_edges']

        # Return a dictionary of combined arrays and their corresponding edge information
        self.combined_data_histogram = {
                'array': combined_array1,
                'x_edges': xedges,
                'y_edges': yedges
        }

        for key in self.histogram_mc_arrays:

            x_title = self.histogram_mc_arrays[key]["title"].split(' ')[0]
            y_title = self.histogram_mc_arrays[key]["title"].split(' ')[1]
            hist_plane= key.split(' ')[1]
            hits = key.split(' ')[6]
            if selected_plane == False:
                plane = hist_plane
            else:
                plane = selected_plane

            if x_title == selected_x_title and y_title == selected_y_title and (hist_plane == plane) and hits == "ALL":

                combined_array2 += self.histogram_mc_arrays[key]["histogram"]
                xedges = self.histogram_mc_arrays[key]['x_edges']
                yedges = self.histogram_mc_arrays[key]['y_edges']

        self.combined_mc_histogram = {
                'array': combined_array2,
                'x_edges': xedges,
                'y_edges': yedges
        }
        
    def calculate_hist_ratio(self):

        xy_DATA_hist = self.combined_data_histogram['array']
        xy_MC_hist = self.combined_mc_histogram['array']
        
        
        data_array = (xy_DATA_hist ) #/ (np.sum(xy_DATA_hist))
        mc_array = (xy_MC_hist ) #/ (np.sum(xy_MC_hist))
        
        ratio = (data_array ) / (mc_array)
        
        self.ratio = ratio

        
    def adjust_both_bins_for_mc_and_data(self):

        min_data_per_bin = self.min_data
        xbins1 = self.combined_data_histogram['x_edges']
        xbins2 = self.combined_mc_histogram['x_edges']
        ybins1 = self.combined_data_histogram['y_edges']
        ybins2 = self.combined_mc_histogram['y_edges']
        array1 = self.combined_data_histogram['array']
        array2 = self.combined_mc_histogram['array']
        
        # calculate the bin centers
        x_centers = (xbins1[:-1] + xbins1[1:]) / 2
        y_centers = (ybins1[:-1] + ybins1[1:]) / 2


        # Combine the bin edges from both histograms for each axis (probably not necessary for these arrays, but this generalizes)
        combined_xbins = np.unique(np.concatenate([xbins1, xbins2]))
        combined_ybins = np.unique(np.concatenate([ybins1, ybins2]))

        # Create a function to determine new bins
        def create_new_bins(combined_bins, min_data_per_bin, data1,data2):
            new_bins = [combined_bins[0]]
            i = 0
            size1 = 0
            size2 = 0
            for edge in combined_bins[1:]:
                # Check if adding this edge would meet the minimum data criterion
                bin_data1 = data1[(data1 >= new_bins[-1]) & (data1 < edge)]
                bin_data2 = data2[(data2 >= new_bins[-1]) & (data2 < edge)]
                size1 += bin_data1.size
                size2 += bin_data2.size
                if (bin_data1.size >= min_data_per_bin and bin_data2.size >=min_data_per_bin):
                    new_bins.append(edge)
                    size1 = 0
                    size2 = 0
            new_bins = np.array(new_bins)
            #if new_bins[-1] != combined_bins[-1]:
               # new_bins[-1] = combined_bins[-1]
            return new_bins
        


        data1_xn = np.sum(array1, axis = 1)
        data1_yn = np.sum(array1, axis = 0)
        data2_xn = np.sum(array2, axis = 1)
        data2_yn = np.sum(array2, axis = 0)
        
        data1_x = np.repeat(x_centers, data1_xn.astype(int))
        data2_x = np.repeat(x_centers, data2_xn.astype(int))
        data1_y = np.repeat(y_centers, data1_yn.astype(int))
        data2_y = np.repeat(y_centers, data2_yn.astype(int))




        # Calculate new bins for each axis
        new_xbins = create_new_bins(combined_xbins, min_data_per_bin, data1_x,data2_x)
        new_ybins = create_new_bins(combined_ybins, min_data_per_bin,  data1_y,data2_y)

        
        new_hist1 = cython_function.redistribute_counts(array1, xbins1, ybins1,new_xbins,new_ybins)
        new_hist2 = cython_function.redistribute_counts(array2, xbins2, ybins2,new_xbins,new_ybins)
        
        def remove_zeros(hist1,hist2,old_xbins,old_ybins):
            xbins = old_xbins
            xzeros = []
            yzeros = []
            ybins = old_ybins
            for i in range(len(xbins) - 1):
                #if len(xbins)  >= 5:
                    row_zeros1 = hist1[i, hist1[i,:] == 0] 
                    row_zeros2 = hist2[i, hist2[i,:] == 0]
                    row_zeros1[:] +=1
                    row_zeros2[:] +=1
                    row_zerosN = row_zeros1.sum() + row_zeros2.sum()
                    if (row_zerosN > 0):
                        xzeros.append(i)
                    break
            
            for j in range(len(ybins) - 1):
                #if len(ybins) >=5:
                    col_zeros1 = hist1[hist1[:,j] == 0, j] 
                    col_zeros2 = hist2[hist2[:,j] == 0, j]
                    col_zeros1[:] +=1
                    col_zeros2[:] +=1
                    col_zerosN = col_zeros1.sum() + col_zeros2.sum()
                    if (col_zerosN > 0):
                        yzeros.append(j)
                        break
                
            xbins = np.delete(xbins,np.array(xzeros, dtype = int))
            ybins = np.delete(ybins,np.array(yzeros, dtype = int))
            hist1 = cython_function.redistribute_counts(hist1, old_xbins, old_ybins, xbins,ybins)
            hist2 = cython_function.redistribute_counts(hist2, old_xbins, old_ybins, xbins,ybins)
            return hist1, hist2, xbins, ybins
                    
                    
        #new_hist1, new_hist2, new_xbins, new_ybins = remove_zeros(new_hist1, new_hist2, new_xbins, new_ybins)
        
        zeros1 = new_hist1[new_hist1 == 0] 
        zeros2 = new_hist2[new_hist2 == 0]
        zeros1[:] +=1
        zeros2[:] +=1
        zerosN = zeros1.sum() + zeros2.sum()  
        while(zerosN > 0):
            new_hist1, new_hist2, new_xbins, new_ybins = remove_zeros(new_hist1, new_hist2, new_xbins, new_ybins)
            zeros1 = new_hist1[new_hist1 == 0] 
            zeros2 = new_hist2[new_hist2 == 0]
            zeros1[:] +=1
            zeros2[:] +=1
            zerosN = zeros1.sum() + zeros2.sum()  
        
        
        self.combined_data_histogram['array'] = new_hist1
        self.combined_mc_histogram['array'] = new_hist2
        self.combined_data_histogram['x_edges'] = new_xbins
        self.combined_mc_histogram['x_edges'] = new_xbins
        self.combined_data_histogram['y_edges'] = new_ybins
        self.combined_mc_histogram['y_edges'] = new_ybins
        
        self.ybins = new_ybins
        self.xbins = new_xbins
        '''
        neg_xbins = self.xbins[self.xbins <= 0 ]
        pos_xbins = self.xbins[self.xbins > 0 ]
        
        
        root_xbins = array.array('d',self.xbins)
        root_pos_xbins = array.array('d',pos_xbins)
        root_neg_xbins = array.array('d',neg_xbins)
        root_ybins = array.array('d',self.ybins)

        self.root_xbins = root_xbins
        self.root_ybins = root_ybins


        xydata1 = ROOT.TH2D(self.X + self.Y + "Positive Data", self.X + self.Y + "Positive Data", nbinsx = int(len(root_pos_xbins)), xbins = root_pos_xbins, nbinsy = int(len(root_ybins)), ybins = root_ybins)

        for i in np.arange(len(root_pos_xbins) - 1):
            for j in np.arange(len(root_ybins) - 1):
                xydata1.Fill(root_pos_xbins[i], root_ybins[j], new_hist1[i,j])
                
        xydata2 = ROOT.TH2D(self.X + self.Y + "Negative Data", self.X + self.Y + "Negative Data", nbinsx = int(len(root_neg_xbins)), xbins = root_neg_xbins, nbinsy = int(len(root_ybins)), ybins = root_ybins)

        for i in np.arange(len(root_neg_xbins) - 1):
            for j in np.arange(len(root_ybins) - 1):
                xydata2.Fill(float(root_neg_xbins[i]), root_ybins[j], new_hist1[i,j])
               
                

            
        self.data_root = xydata1.Add(xydata2)
        
        
        xysim = ROOT.TH2F(self.X + self.Y + " MC", self.X + self.Y + " MC", nbinsx = len(root_xbins), xbins = root_xbins, nbinsy = len(root_ybins), ybins = root_ybins)
        
        for i in np.arange(len(root_xbins) - 1):
            for j in np.arange(len(root_ybins) - 1):
                xysim.Fill(root_xbins[i], root_ybins[j], new_hist2[i,j])
           
        
                
        self.mc_root = xysim
        
        '''
    

    def calculate_average_y_and_std(self,hist):
        if hist == 'mc':
            hist_array = self.combined_mc_histogram['array']
            x_edges = self.combined_mc_histogram['x_edges']
            y_edges = self.combined_mc_histogram['y_edges']
            
            average = np.mean(hist_array, axis=1)
            std = np.std(hist_array, axis=1)
            
            self.mc_y_avg =  average
            self.mc_y_std = std
        elif hist == 'data':
            hist_array = self.combined_data_histogram['array']
            x_edges = self.combined_data_histogram['x_edges']
            y_edges = self.combined_data_histogram['y_edges']
            
            average = np.mean(hist_array, axis=1)
            std = np.std(hist_array, axis=1)
            
            self.data_y_avg =  average
            self.data_y_std = std
        elif hist == 'ratio':
            self.calculate_hist_ratio()
            hist_array = self.ratio
            x_edges = self.combined_data_histogram['x_edges']
            y_edges = self.combined_data_histogram['y_edges']
            
            average = np.mean(hist_array, axis=1)
            std = np.std(hist_array, axis=1)
            
            self.ratio_y_avg =  average
            self.ratio_y_std = std
            
            self.root_x_data_projection = self.data_root.ProfileX()
            self.root_x_mc_projection = self.mc_root.ProfileX()
            
            self.root_ratio = self.root_x_data_projection.Divide(self.root_x_mc_projection)
            
            
            
            
class Histogram_Plotter(Combined_Histograms):
    def __init__(self, selected_x_title, selected_y_title, selected_plane= "All", min_data = 100):
        
        #Initilize class and set up title variables
        
        Combined_Histograms.__init__(self,selected_x_title, selected_y_title, min_data=min_data, selected_plane= selected_plane)
       
        self.X = selected_x_title
        self.Y = selected_y_title
        
        if selected_plane == "All":
            self.plane = 'All Planes'
        else:
            self.plane = "Plane: " + str(selected_plane)
            
        if self.X == ('Theta' or 'Phi'):
            self.xtitle = self.X + " radians"
        else:
            self.xtitle = self.X + " (cm)"
        self.ytitle = self.Y
          
    def plot_stair_avg_y_for_mc_and_data(self):
        self.adjust_x_bins_for_mc_and_data()
        
        self.calculate_average_y_and_std('mc')
        self.calculate_average_y_and_std('data')


        plt.figure(figsize=(10, 6))
        plt.stairs(self.data_y_avg, edges = self.combined_data_histogram['x_edges'], label = 'Data')
        plt.stairs(self.mc_y_avg, edges = self.combined_mc_histogram['x_edges'], label = 'Simulation')
        plt.xlabel(self.xtitle)
        plt.ylabel(self.ytitle)
        plt.title(self.X + " vs. Average " + self.Y)
        plt.legend()
        plt.show()
        
    def plot_2D_hist(self,hist):
        #self.adjust_both_bins_for_mc_and_data()
        
        if hist == 'mc':
            xedges = self.combined_mc_histogram['x_edges']
            yedges = self.combined_mc_histogram['y_edges']
            array = self.combined_mc_histogram['array']
            
            title = self.xtitle + " vs. " + self.ytitle + " Simulation " + self.plane
            
        elif hist == 'data':
            array = self.combined_data_histogram['array']
            xedges = self.combined_data_histogram['x_edges']
            yedges = self.combined_data_histogram['y_edges']
            
            title = self.xtitle + " vs. " + self.ytitle + " Data " + self.plane
            
        elif hist == 'ratio':
            self.calculate_hist_ratio()
            
            array = self.ratio
            xedges = self.combined_data_histogram['x_edges']
            yedges = self.combined_data_histogram['y_edges']
            
            title = self.xtitle + " vs. " + self.ytitle + " Ratios (Data/Simulation) " + self.plane
            

        plt.figure(figsize=(10, 6))
        plt.pcolormesh(xedges, yedges, np.transpose(array), shading = 'flat')
        plt.ylabel(self.ytitle)
        plt.xlabel(self.xtitle)
        plt.title(title)
        plt.colorbar()
        plt.show()

    def plot_average_ratios_for_all_planes(self):
        hist0 = Combined_Histograms(self.X,self.Y,'0',self.min_data)
        hist1 = Combined_Histograms(self.X,self.Y,'1',self.min_data)
        hist2 = Combined_Histograms(self.X,self.Y,'2',self.min_data)

        hist0.adjust_both_bins_for_mc_and_data()
        hist1.adjust_both_bins_for_mc_and_data()
        hist2.adjust_both_bins_for_mc_and_data()
        
        hist0.calculate_average_y_and_std('ratio')
        hist1.calculate_average_y_and_std('ratio')
        hist2.calculate_average_y_and_std('ratio')
        
        x_bin_centers0 = 0.5 * (hist0.combined_data_histogram['x_edges'][:-1] + hist0.combined_data_histogram['x_edges'][1:])
        x_bin_centers1 = 0.5 * (hist1.combined_data_histogram['x_edges'][:-1] + hist1.combined_data_histogram['x_edges'][1:])
        x_bin_centers2 = 0.5 * (hist2.combined_data_histogram['x_edges'][:-1] + hist2.combined_data_histogram['x_edges'][1:])

        
        xerr0 = hist0.combined_data_histogram['x_edges'][1:] - x_bin_centers0
        xerr1 = hist1.combined_data_histogram['x_edges'][1:] - x_bin_centers1
        xerr2 = hist2.combined_data_histogram['x_edges'][1:] - x_bin_centers2 
        
        if self.X == 'Theta':
            mid0 = len(x_bin_centers0) // 2
            mid1 = len(x_bin_centers1) // 2
            mid2 = len(x_bin_centers2) // 2
            
            x_bin_centers0 = x_bin_centers0[:mid0]
            x_bin_centers1 =  x_bin_centers1[:mid1]
            x_bin_centers2 = x_bin_centers2[:mid2]
            
            xerr0 = xerr0[:mid0]
            xerr1 = xerr1[:mid1]
            xerr2 = xerr2[:mid2]
            
            hist0.ratio_y_avg = hist0.ratio_y_avg[:mid0]
            hist1.ratio_y_avg = hist1.ratio_y_avg[:mid1]
            hist2.ratio_y_avg = hist2.ratio_y_avg[:mid2]
            
            hist0.ratio_y_std = hist0.ratio_y_std[:mid0]
            hist1.ratio_y_std = hist1.ratio_y_std[:mid1]
            hist2.ratio_y_std = hist2.ratio_y_std[:mid2]
        
        canvas = ROOT.TCanvas()
        hist0.root_ratio.Draw()
        hist1.root_ratio.Draw()
        hist2.root_ratio.Draw()
        canvas.Update()

        '''
        plt.figure(figsize=(10, 6))
        plt.errorbar(x_bin_centers0 , hist0.ratio_y_avg, yerr=hist0.ratio_y_std, xerr = xerr0 , fmt='o', label='Plane 0')
        plt.errorbar(x_bin_centers1 , hist1.ratio_y_avg, yerr=hist1.ratio_y_std, xerr = xerr1 , fmt='o', label='Plane 1')
        plt.errorbar(x_bin_centers2 , hist2.ratio_y_avg, yerr=hist2.ratio_y_std, xerr = xerr2 , fmt='o', label='Plane 2')
        plt.xlabel(self.xtitle)
        ytitle = self.Y + " Ratio Averages"
        plt.ylabel(ytitle)
        title = self.X + " " + self.Y + " Ratio Averages (Data/Simulation)"
        plt.title(title)
        plt.legend()
        plt.show()
        '''
 
        
class Combined_3D_Histograms(Combined_Histograms):
    def __init__(self, selected_x_title, selected_y_title,selected_z_title, selected_plane= "All", min_data = 100):
        self.X = selected_x_title
        self.Y = selected_y_title
        self.Z = selected_z_title
        self.plane = selected_plane
        if selected_plane == "All":
            self.planetitle = 'All Planes'
        else:
            self.planetitle = selected_plane
        self.min_data = min_data
        
        xzinfo= Combined_Histograms(selected_x_title,selected_z_title, selected_plane, min_data)
        yzinfo  = Combined_Histograms(selected_y_title,selected_z_title, selected_plane, min_data)
        
        self.xbins = xzinfo.combined_data_histogram['x_edges']
        self.zbins = xzinfo.combined_data_histogram['y_edges']
        self.ybins = yzinfo.combined_data_histogram['x_edges']
        
        xzhist = xzinfo.combined_data_histogram['array']
        yzhist = yzinfo.combined_data_histogram['array']
        
        xhist = np.sum(xzhist, axis=1)
        yhist = np.sum(yzhist, axis=1)
        zhist = np.sum(xzhist,axis = 0) + np.sum(yzhist,axis = 0)
        
        xyzdata = ROOT.TH3F(self.X + self.Y + self.Z, self.X + self.Y + self.Z, len(self.xbins), self.xbins[0], self.xbins[-1], len(self.ybins), self.ybins[0], self.ybins[-1], len(self.zbins), self.zbins[0], self.zbins[-1])
        
        for i in np.arange(len(self.xbins) - 1):
            for j in np.arange(len(self.ybins) - 1):
                for k in np.arange(len(self.zbins) - 1):
                    xyzdata.Fill(self.xbins[i], self.ybins[j], self.zbins[k], xhist[i] + yhist[j] + zhist[k])
            
        self.data_3D = xyzdata
        
        xzhist = xzinfo.combined_mc_histogram['array']
        yzhist = yzinfo.combined_mc_histogram['array']
        
        xhist = np.sum(xzhist, axis=1)
        yhist = np.sum(yzhist, axis=1)
        zhist = np.sum(xzhist,axis = 0) + np.sum(yzhist,axis = 0)
        
        xyzmc = ROOT.TH3F(self.X + self.Y + self.Z, self.X + self.Y + self.Z, len(self.xbins), self.xbins[0], self.xbins[-1], len(self.ybins), self.ybins[0], self.ybins[-1], len(self.zbins), self.zbins[0], self.zbins[-1])
        
        for i in np.arange(len(self.xbins) - 1):
            for j in np.arange(len(self.ybins) - 1):
                for k in np.arange(len(self.zbins) - 1):
                    xyzmc.Fill(self.xbins[i], self.ybins[j], self.zbins[k], xhist[i] + yhist[j] + zhist[k])
                    
        self.mc_3D = xyzmc
        
    def projection_2D(self):
        
        #Use Root TProfile fucntion to set up 2d projections of 3d histogram and then store them as numpy arrays
        
        self.root_xy_data_projection = self.data_3D.Project3DProfile( option = "xy")
        
        xedgesN = self.root_xy_data_projection.GetNbinsX()
        yedgesN = self.root_xy_data_projection.GetNbinsY()
        
        xedges = np.array([self.root_xy_data_projection.GetXaxis().GetBinLowEdge(i) for i in range(1, xedgesN + 2)])
        yedges = np.array([self.root_xy_data_projection.GetYaxis().GetBinLowEdge(j) for j in range(1, yedgesN + 2)]) 
        
        self.xy_data_projection = np.array([[self.root_xy_data_projection.GetBinContent(i, j) for j in range(1, xedgesN + 1)] for i in range(1, yedgesN + 1)])[:-1,:-1]
        
        
        self.root_xy_mc_projection = self.mc_3D.Project3DProfile( option = "xy")
        
        xedgesN = self.root_xy_mc_projection.GetNbinsX()
        yedgesN = self.root_xy_mc_projection.GetNbinsY()
        
        xedges = np.array([self.root_xy_mc_projection.GetXaxis().GetBinLowEdge(i) for i in range(1, xedgesN + 2)])
        yedges = np.array([self.root_xy_mc_projection.GetYaxis().GetBinLowEdge(j) for j in range(1, yedgesN + 2)]) 
        
        self.xy_mc_projection = np.array([[self.root_xy_mc_projection.GetBinContent(i, j) for j in range(1, xedgesN + 1)] for i in range(1, yedgesN + 1)])[:-1,:-1]
        
    def plot_projection_ratios(self):
        hists = Combined_Histograms(self.X,self.Y,self.plane,self.min_data)
        
        self.projection_2D()
        
        hists.combined_data_histogram['array'] = self.xy_data_projection
        hists.combined_mc_histogram['array'] = self.xy_mc_projection
        hists.combined_data_histogram['x_edges'] = self.xbins
        hists.combined_data_histogram['y_edges'] = self.ybins
        hists.combined_mc_histogram['x_edges'] = self.xbins
        hists.combined_mc_histogram['y_edges'] = self.ybins
        
        hists.adjust_both_bins_for_mc_and_data()
        
        hists.calculate_hist_ratio()
        
        plt.figure(figsize=(10, 6))
        plt.pcolormesh(hists.combined_data_histogram['x_edges'], hists.combined_data_histogram['y_edges'], hists.ratio.T, shading = 'flat')
        if (self.X in ['Theta', 'Phi']):
            self.xtitle = self.X + " Radians"
        else:
            self.xtitle = self.X + " (cm)"
        if (self.Y in ['Theta' or 'Phi']):
            self.ytitle = self.Y + " Radians"
        else:
            self.ytitle = self.Y + " (cm)"
        plt.ylabel(self.ytitle)
        plt.xlabel(self.xtitle)
        title =  self.X + " vs " + self.Y + " " + self.Z + " Ratios (Data/Simulation) Plane: " + self.planetitle
        plt.title(title)
        plt.colorbar()
        plt.show()
        
        
        
        
        
        