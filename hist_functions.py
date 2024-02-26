import ROOT
import uproot
import numpy as np
import matplotlib.pyplot as plt
import os
import awkward 
from scipy.optimize import curve_fit

def hist_ratio(hists):
    small_constant=1e-6
    xedges = hists['x_edges']
    yedges = hists['y_edges']
    xy_DATA_hist = hists['array_1']
    xy_MC_hist = hists['array_2']
    
    data_array = xy_DATA_hist / (np.sum(xy_DATA_hist))
    mc_array = xy_MC_hist / (np.sum(xy_MC_hist))

    ratio = (data_array + small_constant) / (mc_array + small_constant)
    return ratio

def adjust_bins(histogram, min_data):
    hist = histogram['array']
    xedges = histogram['x_edges']
    yedges = histogram['y_edges']

    new_hist = []
    new_xedges = [xedges[0]]
    cumulative_count = 0

    for i in range(len(hist[0]) - 1):
        cumulative_count += np.sum(hist[i,:])
        if cumulative_count >= min_data:
            new_hist.append(hist[i , :])
            new_xedges.append(xedges[i + 1])
            cumulative_count = 0
        else:
            # Add current column to the next one if the loop is not at the last index
            if i < len(hist[0]) - 2:
                hist[i + 1, :] += hist[i, :]

    # Handle the last bin separately if it hasn't been added
    if cumulative_count > 0:
        new_hist.append(hist[-1,:])
        new_xedges.append(xedges[-1])

    # Convert new_hist to a numpy array
    new_hist = np.array(new_hist)
    new_xedges = np.array(new_xedges)

    return {'array': new_hist, 'x_edges': new_xedges, 'y_edges': yedges}

def average_y_and_std(hist,x,y):
    hist_array = hist
    x_edges = x
    y_edges = y
    
    average_ratio = np.mean(hist_array, axis=1)
    std = np.std(hist_array, axis=1)
    
    return average_ratio, std 

def adjust_bins_for_2_hists(hist1,hist2, min_data):
    new_hist1 = []
    new_hist2 = []
    xedges = hist1['x_edges']
    yedges = hist1['y_edges']
    hist1 = hist1['array']
    hist2 = hist2['array']
    new_xedges = [xedges[0]]
    varhist1 = hist1
    varhist2 = hist2
    cumulative_count1 = 0
    cumulative_count2 = 0

    # Adjusting x bins
    for i in range(len(hist1[0]) - 1):
        cumulative_count1 += np.sum(hist1[i,:])
        cumulative_count2 += np.sum(hist2[i,:])
        if cumulative_count1 >= min_data or cumulative_count2 >= min_data:
            new_hist1.append(varhist1[i , :])
            new_hist2.append(varhist2[i , :])
            new_xedges.append(xedges[i + 1])
            cumulative_count1 = 0
            cumulative_count2 = 0
        else:
            # Add current column to the next one if the loop is not at the last index
            if i < len(hist1[0]) - 2:
                varhist1[i + 1, :] += varhist1[i, :]
                varhist2[i + 1, :] += varhist2[i, :]

    # Handle the last bin separately if it hasn't been added
    if cumulative_count1 > 0 and cumulative_count2 > 0:
        new_hist1.append(hist1[-1,:])
        new_hist2.append(hist2[-1,:])
        new_xedges.append(xedges[-1])
    
     # Convert new_hist to a numpy array
    new_hist1 = np.array(new_hist1)
    new_hist2 = np.array(new_hist2)
    new_xedges = np.array(new_xedges)

    return {'array_1': new_hist1, 'array_2': new_hist2, 'x_edges': new_xedges, 'y_edges': yedges}
    
def adjust_both_bins_for_2_hists(hist1,hist2, xedges, yedges, min_data):
    new_hist1 = []
    new_hist2 = []
    new_xedges = [xedges[0]]
    new_yedges = [yedges[0]]
    varhist1 = hist1
    varhist2 = hist2
    cumulative_count1 = 0
    cumulative_count2 = 0

    # Adjusting x bins
    for i in range(len(hist1[0]) - 1):
        cumulative_count1 += np.sum(hist1[i,:])
        cumulative_count2 += np.sum(hist2[i,:])
        if cumulative_count1 >= min_data and cumulative_count2 >= min_data:
            new_hist1.append(varhist1[i , :])
            new_hist2.append(varhist2[i , :])
            new_xedges.append(xedges[i + 1])
            cumulative_count1 = 0
            cumulative_count2 = 0
        else:
            # Add current column to the next one if the loop is not at the last index
            if i < len(hist1[0]) - 2:
                varhist1[i + 1, :] += varhist1[i, :]
                varhist2[i + 1, :] += varhist2[i, :]

    # Handle the last bin separately if it hasn't been added
    if cumulative_count1 > 0 or cumulative_count2 > 0:
        new_hist1.append(hist1[-1,:])
        new_hist2.append(hist2[-1,:])
        new_xedges.append(xedges[-1])
    
     # Convert new_hist to a numpy array
    new_hist1 = np.array(new_hist1)
    new_hist2 = np.array(new_hist2)
    
    # Now for y
    cumulative_count1 = 0
    cumulative_count2 = 0
    for i in range(len(new_hist1[1]) - 1):
        if i < len(new_hist1[1]) - 2:
            cumulative_count1 += np.sum(new_hist1[:,i])
            cumulative_count2 += np.sum(new_hist2[:,i])
            if cumulative_count1 >= min_data and cumulative_count2 >= min_data:
                new_yedges.append(yedges[i + 1])
                cumulative_count1 = 0
                cumulative_count2 = 0
            else:
                if i < len(hist1[1]) - 2:
                    new_hist1[:, i + 1] += new_hist1[:, i]
                    new_hist2[:, i + 1] += new_hist2[:, i]

    if cumulative_count1 > 0 or cumulative_count2 > 0:
        new_yedges.append(yedges[-1])
        
    # Convert edges to a numpy array
    new_xedges = np.array(new_xedges)
    new_yedges = np.array(new_yedges)
    
    # y bins that were combined
    ybins = np.setdiff1d(yedges, new_yedges)
    indices = np.where(np.isin(yedges,ybins)) 
    ones = np.ones(shape = np.shape(indices), dtype= int )
    indices = indices - ones
    
    new_hist1 = np.delete(new_hist1, indices, axis = 1)
    new_hist2 = np.delete(new_hist2, indices, axis = 1)
    
    return {'array_1': new_hist1, 'array_2': new_hist2, 'x_edges': new_xedges, 'y_edges': new_yedges}
 
class HistogramAnalyzer:
    def __init__(self, file_data_path, file_mc_path, small_constant=1e-6):
        self.file_data_path = file_data_path
        self.file_mc_path = file_mc_path
        self.small_constant = small_constant
        self.ratio_histograms = {}
        
    def normalize_2d_hist(self, hist, x_edges, y_edges):
        area = (x_edges[-1] - x_edges[0]) * (y_edges[-1] - y_edges[0])
        return hist / (hist.sum() * area)
    
    def calculate_histogram_ratios(self):
        file_data = uproot.open(self.file_data_path)
        file_mc = uproot.open(self.file_mc_path)

        for key in file_data.keys():
            if key in file_mc:
                obj_data = file_data[key]
                obj_mc = file_mc[key]

                if 'TH2' in obj_data.classname and 'TH2' in obj_mc.classname:
                    data_array, data_x_edges, data_y_edges = obj_data.to_numpy()
                    mc_array, mc_x_edges, mc_y_edges = obj_mc.to_numpy()

                    data_array = self.normalize_2d_hist(data_array, data_x_edges, data_y_edges)
                    mc_array = self.normalize_2d_hist(mc_array, mc_x_edges, mc_y_edges)

                    ratio = (data_array + self.small_constant) / (mc_array + self.small_constant)

                    # Retrieve the axes titles from the original histograms
                    x_title = obj_data.title.split(' ')[0]
                    y_title = obj_data.title.split(' ')[1]

                    self.ratio_histograms[key] = {
                    'ratio': ratio,
                    'x_edges': data_x_edges,
                    'y_edges': data_y_edges,
                    'x_title': x_title,
                    'y_title': y_title
                    }
        return self.ratio_histograms

    def save_histograms_to_root_file(self, output_file_name):
        root_file = ROOT.TFile(output_file_name, "RECREATE")

        for hist_name, data in self.ratio_histograms.items():
            ratio = data['ratio']
            x_edges = data['x_edges']
            y_edges = data['y_edges']
            x_title = data['x_title']
            y_title = data['y_title']
            
            x_array = ROOT.std.vector('double')()
            [x_array.push_back(edge) for edge in x_edges]

            y_array = ROOT.std.vector('double')()
            [y_array.push_back(edge) for edge in y_edges]

            hist = ROOT.TH2D(hist_name, hist_name[:-6] + ";" + x_title + ";" + y_title, len(x_edges)-1, x_array.data(), len(y_edges)-1, y_array.data())

            for i in range(len(x_edges)-1):
                for j in range(len(y_edges)-1):
                    hist.SetBinContent(i + 1, j + 1, ratio[i, j])

            hist.Write()

        root_file.Close()

    def plot_and_save_histograms_by_selected_axes(self, selected_x_values, selected_y_values, output_folder, plots_per_figure=9):
        # Ensure the output directory exists
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        # Group histograms based on the selected x and y axis values
        grouped_histograms = {}
        for hist_name, data in self.ratio_histograms.items():
            x_axis_label = data['x_title']
            y_axis_label = data['y_title']

            if x_axis_label in selected_x_values and y_axis_label in selected_y_values:
                axis_label_pair = (x_axis_label, y_axis_label)

                if axis_label_pair not in grouped_histograms:
                    grouped_histograms[axis_label_pair] = []
                grouped_histograms[axis_label_pair].append((hist_name, data))

        # Initialize figure index
        figure_index = 1

        for axis_labels, histograms in grouped_histograms.items():
            for i, (hist_name, data) in enumerate(histograms):
                if i % plots_per_figure == 0:
                    if i > 0:  # Save the previous figure
                        plt.tight_layout()
                        plt.savefig(os.path.join(output_folder, f'Grouped_Histograms_{figure_index}.png'))
                        plt.close()
                        figure_index += 1

                    # Create a new figure
                    fig, axs = plt.subplots(plots_per_figure // 3, 3, figsize=(15, 10))
                    axs = axs.flatten()  # Flatten to 1D array for easy indexing

                # Get current axis
                ax = axs[i % plots_per_figure]

                # Plot the histogram ratio
                ratio = data['ratio']
                x_edges = data['x_edges']
                y_edges = data['y_edges']
                pcm = ax.pcolormesh(x_edges, y_edges, ratio.T, shading='auto', cmap='viridis')
                fig.colorbar(pcm, ax=ax, orientation='vertical')

                # Set titles and labels
                ax.set_title(hist_name)
                ax.set_xlabel(data['x_title'])
                ax.set_ylabel(data['y_title'])

            # Turn off any unused subplots in the last figure
            for j in range(i % plots_per_figure + 1, len(axs)):
                axs[j].axis('off')

            # Save the last figure
            plt.tight_layout()
            plt.savefig(os.path.join(output_folder, f'Grouped_Histograms_{figure_index}.png'))
            plt.close()

    def combine_histograms_by_axes_and_plane(self, file_path, selected_plane, selected_x_title, selected_y_title):
        xedges = np.zeros(shape = (100))
        yedges = np.zeros(shape = (100))
        combined_array = np.zeros(shape = (100,100))
        file = uproot.open(file_path)

        for key in file.keys():
            obj = file[key]

            if 'TH2' in obj.classname:  # Check for 2D histogram
                # Retrieve the histogram data and axis information
                hist_array, x_edges, y_edges = obj.to_numpy()
                x_title = obj.title.split(' ')[0]
                y_title = obj.title.split(' ')[1]
                hist_plane= key.split(' ')[1]
                hits = key.split(' ')[6]

                # Check if the axes titles match the selected ones
                if x_title == selected_x_title and y_title == selected_y_title and (hist_plane == selected_plane) and hits == "ALL":

                    # Add the current histogram to the combined array for this pair of axes
                    combined_array += hist_array
                    xedges = x_edges
                    yedges = y_edges

        # Return a dictionary of combined arrays and their corresponding edge information
        combined_histograms_info = {
                'array': combined_array,
                'x_edges': xedges,
                'y_edges': yedges
        }

        return combined_histograms_info
         
    def plot_average_y_by_title_part(self, combined_histograms, selected_plane, selected_x_title, selected_y_title):
        hist_array = combined_histograms[(selected_x_title, selected_y_title)]['array']
        x_edges = combined_histograms[(selected_x_title, selected_y_title)]['x_edges']
        y_edges = combined_histograms[(selected_x_title, selected_y_title)]['y_edges']

        # Calculate the y-bin centers
        y_bin_centers = 0.5 * (y_edges[1:] + y_edges[:-1])

        # Calculate the sum of weights (counts) for each x-bin
        sum_weights = hist_array.sum(axis=0)  # sum over y

        # Filter out bins with fewer than 10 entries
        valid_bins = sum_weights >= 10

        # Calculate the average y-value for each valid x-bin
        with np.errstate(divide='ignore', invalid='ignore'):
            weighted_y_values = hist_array * y_bin_centers[:, None]  # weight y-values by their counts
            sum_weighted_y = weighted_y_values.sum(axis=0)  # sum weighted y-values for each x-bin
            average_y = np.divide(sum_weighted_y, sum_weights, where=valid_bins)  # compute average y

            # Calculate standard deviation for each valid x-bin
            squared_diffs = (y_bin_centers[:, None] - average_y)**2 * hist_array
            sum_squared_diffs = squared_diffs.sum(axis=0)
            std_dev_y = np.sqrt(np.divide(sum_squared_diffs, sum_weights, where=valid_bins))

        # Calculate the error on the average y for valid bins
        errors_on_mean = np.divide(std_dev_y, np.sqrt(sum_weights), where=valid_bins)

        # Adjust x-bin edges for positive bins
        positive_bin_indices = np.where(x_edges > 0)[0]
        for i in positive_bin_indices[:-1]:  # Skip the last edge
            # Double the distance to the next edge for positive bins
            x_edges[i+1] = x_edges[i] + 2 * (x_edges[i+1] - x_edges[i])

        # Recalculate x-bin centers after modifying edges
        x_bin_centers = 0.5 * (x_edges[:-1] + x_edges[1:]) 

        # Ensure we only use the centers for valid bins
        x_bin_centers = x_bin_centers[valid_bins]
        average_y = average_y[valid_bins]
        errors_on_mean = errors_on_mean[valid_bins]

        plt.figure(figsize=(10, 6))
        plt.errorbar(x_bin_centers , average_y, yerr=errors_on_mean, fmt='o', label=f'{selected_plane}')
        plt.xlabel(selected_x_title)
        plt.ylabel(f'Average {selected_y_title}')
        plt.title(f'Average {selected_y_title} vs {selected_x_title} for {selected_plane}')
        plt.legend()
        plt.grid(True)
        plt.show()

    def combine_histograms_by_axes(self,file_path,selected_x_title,selected_y_title): 
        combined_histograms = np.zeros(shape = (100,100))
        xedges = np.zeros(shape = (100))
        yedges = np.zeros(shape = (100))
        file = uproot.open(file_path)

        for key in file.keys():
            obj = file[key]

            if 'TH2' in obj.classname:  # Check for 2D histogram
                # Retrieve the histogram data and axis information
                hist_array, x_edges, y_edges = obj.to_numpy()
                x_title = obj.title.split(' ')[0]
                y_title = obj.title.split(' ')[1]
                hits = key.split(' ')[6]

                # Check if the axes titles match the selected ones
                if x_title == selected_x_title and y_title == selected_y_title and hits == "ALL":

                    # Add the current histogram to the combined array for this pair of axes
                    combined_histograms += hist_array
                    xedges = x_edges
                    yedges = y_edges
                    

        # Return a dictionary of combined arrays and their corresponding edge information
        combined_arrays_info = {
                'array': combined_histograms,
                'x_edges': xedges,
                'y_edges': yedges
        }
        

        return combined_arrays_info
    
    def figure_5_plot(self,X,Y,bin_min):
        xAEx = self.combine_histograms_by_axes('histogramsDATA.root', X, Y)
        xAMC = self.combine_histograms_by_axes('histogramsMC.root', X, Y)

        xAEx = adjust_bins(xAEx, bin_min)
        xAMC = adjust_bins(xAMC, bin_min)

        AvgAEx = average_y(xAEx['array'],xAEx['x_edges'],xAEx['y_edges'])
        AvgAMC = average_y(xAMC['array'],xAMC['x_edges'],xAMC['y_edges'])

        xExBins = xAEx['x_edges']
        xMCBins = xAMC['x_edges']

        plt.figure(figsize=(10, 6))
        plt.stairs(AvgAEx, edges = xExBins, label = 'Data')
        plt.stairs(AvgAMC, edges = xMCBins, label = 'Simulation')
        plt.xlabel(X)
        plt.ylabel(Y)
        plt.legend()
        plt.show()
        
    def figure_7_plot(self,file,X,Y):
        
        hist = self.combine_histograms_by_axes(file, X, Y)
        
        plt.figure(figsize=(10, 6))
        plt.pcolormesh(hist['x_edges'], hist['y_edges'], np.transpose(hist['array']), shading = 'flat')
        plt.ylabel(Y)
        plt.xlabel(X)
        plt.show()

    def figure_9_plot(self,X,Y,min_count):
        MC_hist0 = self.combine_histograms_by_axes_and_plane('histogramsMC.root','0',X,Y)
        MC_hist1 = self.combine_histograms_by_axes_and_plane('histogramsMC.root','1',X,Y)
        MC_hist2 = self.combine_histograms_by_axes_and_plane('histogramsMC.root','2',X,Y)
        DATA_hist0 = self.combine_histograms_by_axes_and_plane('histogramsDATA.root','0',X,Y)
        DATA_hist1 = self.combine_histograms_by_axes_and_plane('histogramsDATA.root','1',X,Y)
        DATA_hist2 = self.combine_histograms_by_axes_and_plane('histogramsDATA.root','2',X,Y)
        
        #adjust bins for min count
        hist0 = adjust_bins_for_2_hists(DATA_hist0, MC_hist0, min_count)
        hist1 = adjust_bins_for_2_hists(DATA_hist1, MC_hist1, min_count)
        hist2 = adjust_bins_for_2_hists(DATA_hist2, MC_hist2, min_count)
        
        avg_ratio0, std0 = average_y_and_std(hist_ratio(hist0),hist0['x_edges'],hist0['y_edges'])
        avg_ratio1, std1 = average_y_and_std(hist_ratio(hist1),hist1['x_edges'],hist1['y_edges'])
        avg_ratio2, std2 = average_y_and_std(hist_ratio(hist2),hist2['x_edges'],hist2['y_edges'])
        
        x_bin_centers0 = 0.5 * (hist0['x_edges'][:-1] + hist0['x_edges'][1:])
        x_bin_centers1 = 0.5 * (hist1['x_edges'][:-1] + hist1['x_edges'][1:]) 
        x_bin_centers2 = 0.5 * (hist2['x_edges'][:-1] + hist2['x_edges'][1:])
        
        xerr0 = hist0['x_edges'][1:] - x_bin_centers0
        xerr1 = hist1['x_edges'][1:] - x_bin_centers1
        xerr2 = hist2['x_edges'][1:] - x_bin_centers2 
        
        if X == 'Theta':
            mid0 = len(x_bin_centers0) // 2
            mid1 = len(x_bin_centers1) // 2
            mid2 = len(x_bin_centers2) // 2
            
            x_bin_centers0 = x_bin_centers0[:mid0]
            x_bin_centers1 =  x_bin_centers1[:mid1]
            x_bin_centers2 = x_bin_centers2[:mid2]
            
            xerr0 = xerr0[:mid0]
            xerr1 = xerr1[:mid1]
            xerr2 = xerr2[:mid2]
            
            avg_ratio0 = avg_ratio0[:mid0]
            avg_ratio1 = avg_ratio1[:mid1]
            avg_ratio2 = avg_ratio2[:mid2]
            
            std0 = std0[:mid0]
            std1 = std1[:mid1]
            std2 = std2[:mid2]

        
        plt.figure(figsize=(10, 6))
        plt.errorbar(x_bin_centers0 , avg_ratio0, yerr=std0, xerr = xerr0 , fmt='o', label='Plane 0')
        plt.errorbar(x_bin_centers1 , avg_ratio1, yerr=std1, xerr = xerr1 , fmt='o', label='Plane 1')
        plt.errorbar(x_bin_centers2 , avg_ratio2, yerr=std2, xerr = xerr2 , fmt='o', label='Plane 2')
        plt.xlabel(X)
        ytitle = Y + " Ratio Averages"
        plt.ylabel(ytitle)
        plt.legend()
        plt.show()
        
    def create_combined_3D_ROOT_hist(self,file_path, plane, X, Y, Z,bin_min):
        xzhist, xbins, zbins = self.combine_histograms_by_axes_and_plane(file_path,plane, X,Z).values()
        yzhist, ybins, _  = self.combine_histograms_by_axes_and_plane(file_path, plane, Y,Z).values()
        
        xhist = np.sum(xzhist, axis=1)
        yhist = np.sum(yzhist, axis=1)
        zhist = np.sum(xzhist,axis = 0) + np.sum(yzhist,axis = 0)
        
        xyzroot = ROOT.TH3F(X + Y + Z, X + Y + Z, len(xbins), xbins[0], xbins[-1], len(ybins), ybins[0], ybins[-1], len(zbins), zbins[0], zbins[-1])
        
        for i in np.arange(len(xbins) - 1):
            for j in np.arange(len(ybins) - 1):
                for k in np.arange(len(zbins) - 1):
                    xyzroot.Fill(xbins[i], ybins[j], zbins[k], xhist[i] + yhist[j] + zhist[k])
            
        return xyzroot
                
    def average_3D_ROOT_hist(self,xyzroot):
        xyzprofile = xyzroot.Project3DProfile( option = "xy")
        
        return xyzprofile
    
    def plot_ratios_3D_Average(self,plane, X,Y,Z,min_data): 
        
        dataroot = self.average_3D_ROOT_hist(self.create_combined_3D_ROOT_hist('histogramsDATA.root',plane, X, Y, Z,min_data))
        
        xedgesN = dataroot.GetNbinsX()
        yedgesN = dataroot.GetNbinsY()
        
        xedges = np.array([dataroot.GetXaxis().GetBinLowEdge(i) for i in range(1, xedgesN + 2)])
        yedges = np.array([dataroot.GetYaxis().GetBinLowEdge(j) for j in range(1, yedgesN + 2)]) 
        
        dataint = dataroot.Integral("width")
        
        datahist = np.array([[dataroot.GetBinContent(i, j) for j in range(1, xedgesN + 1)] for i in range(1, yedgesN + 1)])
        
        mcroot = self.average_3D_ROOT_hist(self.create_combined_3D_ROOT_hist('histogramsMC.root',plane, X, Y, Z,min_data))
        mchist = np.array([[mcroot.GetBinContent(i, j) for j in range(1, xedgesN + 1)] for i in range(1, yedgesN + 1)])
        mcint = mcroot.Integral("width")
        
        hists = adjust_both_bins_for_2_hists(datahist, mchist, xedges, yedges,min_data)
        
        #hists = {'array_1': datahist, 'array_2': mchist, 'x_edges': xedges, 'y_edges': yedges}
        
        ratio = hist_ratio(hists)
        
        #ratio = (datahist/np.sum(datahist) + self.small_constant) /( mchist/np.sum(mchist) + self.small_constant)
        
        plt.figure(figsize=(10, 6))
        plt.pcolormesh(hists['y_edges'],hists['x_edges'], ratio, shading = 'flat')
        plt.xlabel(X + " (cm)")
        plt.ylabel(Y+ " (cm)")
        name = Z + ' Ratios (Data/Simulation) for Plane ' + plane 
        plt.title(name)
        plt.colorbar()
        plt.show()

        
        
        
        
        
 

