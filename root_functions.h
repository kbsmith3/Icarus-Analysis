#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <vector>
#include <string>
#include <TH2.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TH1.h>
#include <TH3D.h>
#include <TH3.h>
#include <TKey.h>
#include <TCanvas.h>
#include <TApplication.h> 
#include <TProfile.h>
#include <TAxis.h>
#include <TArrayD.h>
#include <algorithm>
#include <memory> 
#include <TProfile2D.h>
#include <random>
#include <stdexcept>
#include <unordered_map>

using namespace std;

//=============================================================================================================================

//This function loops through root file and returns the ALL Hits version every histgoram with a given x and y axis
vector<unique_ptr<TH2D>> loop_through_root_file(TFile* file, const string& X, const string& Y) {
    vector<unique_ptr<TH2D>> histograms;
    if (!file || file->IsZombie()) {
        cerr << "Failed to open file " << endl;
        return histograms;
    }
    TList* list = file->GetListOfKeys();
    if (!list) {
        cerr << "No keys found in file" << endl;
        file->Close();
        delete file;
        return histograms;
    }
    TIter iter(list);
    TKey* key;
    while ((key = (TKey*)iter())) {
        string keyName = key->GetName();
        const string& hits = "ALL";
        if (keyName.find(X) != string::npos && keyName.find(Y) != string::npos && keyName.find(hits) != string::npos) {
            TObject* obj = key->ReadObj();
            if (obj && obj->IsA()->InheritsFrom("TH2D")) { // Ensure the object is a TH2D
                TH2D* clonedHist = static_cast<TH2D*>(obj->Clone());
                histograms.push_back(unique_ptr<TH2D>(clonedHist)); // Add the cloned histogram to the vector
            }
        }
    }  
    return histograms;
}

//=============================================================================================================================

//Combines all histograms with certain string in name; used to combine the planes
void combineHistogramsWithName(vector<unique_ptr<TH2D>>& histograms, const string& nameContains, unique_ptr<TH2D>& combinedHistogram) {
    for (const auto& hist : histograms) {
        string thisname = hist->GetName();
        if (thisname.find(nameContains) != string::npos) {
            if (!combinedHistogram) {
                // For the first matching histogram, clone it to initialize the combinedHistogram
                combinedHistogram.reset((TH2D*)hist->Clone((string(hist->GetName()) + "_combined").c_str()));
            } else {
                // Add the current histogram to the combinedHistogram
                combinedHistogram->Add(hist.get());
            }
        }
    }
    if (combinedHistogram) {
        //cout << "Combined histogram: " << combinedHistogram->GetName() << endl;
    } else {
        cout << "No histograms matching '" << nameContains << "' were found to combine." << endl;
    }
}

//=============================================================================================================================

//creates new bins based on the minimum information needed 
TArrayD* create_new_bins(TArrayD* bins, TH1D* data_prof, TH1D* mc_prof, const double& min_data) {
    vector<double> new_bins;
    double data_N = 0;
    double mc_N = 0;

    //add first bin
    new_bins.push_back((*bins)[0]);

    // only add bins if there is enough data
    for (Int_t i = 0; i < bins->GetSize(); ++i) {
        data_N += data_prof->GetBinContent(i + 1);
        mc_N += mc_prof->GetBinContent(i + 1);
        if ((data_N >= min_data) && (mc_N >= min_data)) {
            new_bins.push_back((*bins)[i]); 
            data_N = 0;
            mc_N = 0;
        }
    }

    TArrayD* bin_array = new TArrayD(new_bins.size()); 
    for (size_t i = 0; i < new_bins.size(); ++i) {
        bin_array->AddAt(new_bins[i], i);
    }

    return bin_array;
}

//=============================================================================================================================

// Function to redistribute counts from an old 2D histogram to a new binning scheme using TArrayD
TH2D* RedistributeCounts(TH2D* old_hist, TArrayD* x_edges_new, TArrayD* y_edges_new) {
    // Create a new histogram with the new binning scheme
    TH2D* new_hist =  new TH2D("new_hist", "Redistributed Histogram",
                                           x_edges_new->GetSize() - 1, x_edges_new->GetArray(),
                                           y_edges_new->GetSize() - 1, y_edges_new->GetArray());

    // Iterate through all old bins
    for (int i = 1; i <= old_hist->GetNbinsX(); ++i) {
        for (int j = 1; j <= old_hist->GetNbinsY(); ++j) {
            double x_center = old_hist->GetXaxis()->GetBinCenter(i);
            double y_center = old_hist->GetYaxis()->GetBinCenter(j);
            double content = old_hist->GetBinContent(i, j);

            // Assuming the old bins map directly to new bins or fall entirely within new bins
            // Find the new bin indices for the given bin center
            int new_i = new_hist->GetXaxis()->FindBin(x_center);
            int new_j = new_hist->GetYaxis()->FindBin(y_center);

            // Check if the bin center is within the bounds of the new histogram
            if (new_i > 0 && new_i <= new_hist->GetNbinsX() && new_j > 0 && new_j <= new_hist->GetNbinsY()) {
                // Add the old bin's content to the corresponding new bin
                double existingContent = new_hist->GetBinContent(new_i, new_j);
                new_hist->SetBinContent(new_i, new_j, existingContent + content);
            }
        }
    }
    return new_hist;
}

//=============================================================================================================================

//For creating the bins of a histogram that was inintiallized with non variable bin widths 
TArrayD* CreateUniformBinArray(int nBins, double xMin, double xMax) {
    TArrayD* bins = new TArrayD(nBins + 1); // nBins + 1 because there are nBins + 1 edges for nBins
    double binWidth = (xMax - xMin) / nBins;
    for (int i = 0; i <= nBins; ++i) {
        bins->AddAt(xMin + i * binWidth,i);
    }
    return bins;
}

//=============================================================================================================================

//finds the mode of a vector or picks a random element tied with highest occurances if there is not one mode
template<typename T>
T findMode(const std::vector<T>& vec) {
    if (vec.empty()) throw std::runtime_error("Empty vector has no mode");

    std::unordered_map<T, int> countMap;
    // Count the occurrences of each element
    for (const auto& elem : vec) {
        ++countMap[elem];
    }

    // Find the maximum frequency
    int maxFrequency = std::max_element(countMap.begin(), countMap.end(),
                                        [](const std::pair<T, int>& a, const std::pair<T, int>& b) {
                                            return a.second < b.second;
                                        })->second;

    // Collect all elements that are tied for the maximum frequency
    std::vector<T> modes;
    for (const auto& pair : countMap) {
        if (pair.second == maxFrequency) {
            modes.push_back(pair.first);
        }
    }

    // Seed with a real random value, if available
    std::random_device rd;

    // Choose a random index
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib(0, modes.size() - 1);
    int randomIndex = distrib(gen);

    // Return a random element from those tied with the maximum frequency
    return modes[randomIndex];
}

//=============================================================================================================================

//checks if there are any bins with zero entries
bool hasZeroBinEntries(TH2D* hist) {
    for (int i = 1; i <= hist->GetNbinsX() ; ++i) {
        for (int j = 1; j <= hist->GetNbinsY() ; ++j) {
            if (hist->GetBinContent(i, j) == 0) {
                return true; // Found a bin with zero entries
            }
        }
    }
    return false; // No bins with zero entries
}

//=============================================================================================================================

//same as above but for 3D hist
bool hasZeroBinEntries3D(TH3D* hist) {
    for (int i = 1; i <= hist->GetNbinsX() ; ++i) {
        for (int j = 1; j <= hist->GetNbinsY() ; ++j) {
            for (int k = 1; k <= hist->GetNbinsZ() ; ++k){
                if (hist->GetBinContent(i, j, k) == 0) {
                    return true; // Found a bin with zero entries
                }
            }
        }
    }
    return false; // No bins with zero entries
}

//=============================================================================================================================

//function for taking a histogram and adjusting the bins so that there are no entries with zero
vector<TArrayD*> remove_zeros( TH2D* old_hist1, TH2D* old_hist2) {
    // Access bin arrays
    int nxbins = old_hist1->GetXaxis()->GetNbins();
    int nybins = old_hist1->GetYaxis()->GetNbins();

    //make vector of bins because they are easy to work with
    vector<double> xbins, ybins;

    // Fill the vectors with the bin edges
    for (int i = 0; i <= nxbins - 1; i++) {
        xbins.push_back(old_hist1->GetXaxis()->GetBinLowEdge(i+1));
    }

    for (int j = 0; j <= nybins - 1; j++) {
        ybins.push_back(old_hist1->GetYaxis()->GetBinLowEdge(j+1));
    }

    // Clone the original histogram for manipulation
    TH2D* new_hist1 = new TH2D(*old_hist1);
    TH2D* new_hist2 = new TH2D(*old_hist2);

    bool zeros1 = hasZeroBinEntries(old_hist1);
    bool zeros2 = hasZeroBinEntries(old_hist2);

    // Create vector that will contain everytime there is a bin with zero in it, then the x and y bins with the most amount of zeros will be removed 
    vector<int> xzeros, yzeros;

    //add count to stop infinte loop
    int count = 0;

    //loops thorugh histogram and removes most the x and y bins with the most zeros until no zeros remain (or until nxbins is reached because this loop can go on forever and that seems like a good cutoff)
    while (zeros1 || zeros2){
        count += 1;
        for (int i = 0 ; i < xbins.size(); i++){ 
            for (int j = 0; j <  ybins.size(); j++){
                Double_t content1 = new_hist1->GetBinContent(i + 1, j + 1);
                Double_t content2 = new_hist2->GetBinContent(i + 1, j + 1);
                if (content1 == 0 || content2 == 0){
                    xzeros.push_back(i);
                    yzeros.push_back(j);
                }
            }
        }
        int i = findMode(xzeros);
        int j = findMode(yzeros);
        //remove bins with the most zeros 
        if (xbins.size() > 5){
            xbins.erase(xbins.begin() + i - 1);
        }
        if (ybins.size() > 5){
            ybins.erase(ybins.begin()+ j - 1);
        }
         // Creating a new histogram to renter the loop with the bins removed
        //make a new array 
        TArrayD* xbin_array = new TArrayD(xbins.size(), xbins.data());
        TArrayD* ybin_array = new TArrayD(ybins.size(), ybins.data());

        new_hist1 = RedistributeCounts(new_hist1, xbin_array, ybin_array);
        new_hist2 = RedistributeCounts(new_hist2, xbin_array, ybin_array);

        delete xbin_array;
        delete ybin_array;

        xzeros.clear();
        yzeros.clear();

        zeros1 = hasZeroBinEntries(new_hist1);
        zeros2 = hasZeroBinEntries(new_hist2);

        //stop code from running forever if zeros are not removed
        if (count > nxbins){
            cout<< "Can't Remove all Zeros" << endl;
            break;
        }
    }
    //return bins 
    TArrayD* xbin_array = new TArrayD(xbins.size(), xbins.data());
    TArrayD* ybin_array = new TArrayD(ybins.size(), ybins.data());

    vector<TArrayD*> bins = {xbin_array,ybin_array};
    return bins;
}

//=============================================================================================================================

//finds ratio of TProfiles by Divide but manuelly calculating errors
TProfile* ratio_profile(TProfile* profileA, TProfile* profileB){
    auto ratioProfile = new TProfile(*profileA); // Clone profileA to keep binning the same

    //set contents to ratio 
    ratioProfile->Divide(profileB);

    int nBins = profileA->GetNbinsX(); // Assuming both profiles have the same number of bins

    for (int iBin = 1; iBin <= nBins; ++iBin) {
        double A = profileA->GetBinContent(iBin);
        double B = profileB->GetBinContent(iBin);
        double sigmaA = profileA->GetBinError(iBin);
        double sigmaB = profileB->GetBinError(iBin);
        if (B != 0) { // Avoid division by zero
            double R = A / B;
            double sigmaR = R * sqrt(pow(sigmaA / A, 2) + pow(sigmaB / B, 2));
            
            ratioProfile->SetBinError(iBin, sigmaR);
        } else {
            // Handle division by zero if necessary, e.g., by setting content to 0 and error to a large number
            ratioProfile->SetBinContent(iBin, 0);
            ratioProfile->SetBinError(iBin, std::numeric_limits<double>::max());
        }
    }
    return ratioProfile;
}

//=============================================================================================================================

  // ----------------------------------------------------------------
 //Creates a new vector of histograms with bins adjusted for minimum data and no zero entries  while making sure both histograms have the same bins({Data_hist, MC_hist})
// Uses many of the above functions
vector<TH2D*> adjust_bins( unique_ptr<TH2D>& data_hist, unique_ptr<TH2D>& mc_hist, const double& min_data){
     //----------------------------------------------------------------
    // create bins for histograms that don't have array of them
    int nBinsX = data_hist->GetXaxis()->GetNbins();
    double xMin = data_hist->GetXaxis()->GetXmin();
    double xMax = data_hist->GetXaxis()->GetXmax();

    int nBinsY = data_hist->GetYaxis()->GetNbins();
    double yMin = data_hist->GetYaxis()->GetXmin();
    double yMax = data_hist->GetYaxis()->GetXmax();

    // Create the TArrayD for the x-axis bins
    TArrayD* xbins = CreateUniformBinArray(nBinsX, xMin, xMax);
    TArrayD* ybins = CreateUniformBinArray(nBinsY, yMin, yMax);
     //----------------------------------------------------------------
    // adjusts the x and y bins for the minimum amount of data needed
    TH1D* data_x_proj = data_hist->ProjectionX();
    TH1D* data_y_proj = data_hist->ProjectionY();
    TH1D* mc_x_proj = mc_hist->ProjectionX();
    TH1D* mc_y_proj = mc_hist->ProjectionY();

    auto new_xbins = create_new_bins(xbins,data_x_proj,mc_x_proj,min_data);
    auto new_ybins = create_new_bins(ybins,data_y_proj,mc_y_proj,min_data);

    //delete projections
    TH2D* data_hist_var = data_hist.get();
    TH2D* mc_hist_var = mc_hist.get();

    // Use the redistributed counts to create new histograms
    auto new_data_hist = RedistributeCounts(data_hist_var, new_xbins, new_ybins);
    auto new_mc_hist = RedistributeCounts(mc_hist_var, new_xbins, new_ybins);
    // ----------------------------------------------------------------
    //remove the zeros from the histogram
    //auto newest_bins = remove_zeros(new_data_hist,new_mc_hist);
    
    //new_data_hist = RedistributeCounts(new_data_hist, newest_bins[0], newest_bins[1]);
    //new_mc_hist = RedistributeCounts(new_mc_hist, newest_bins[0], newest_bins[1]);

    // ----------------------------------------------------------------
    // Create a vector to hold the new histograms
    std::vector<TH2D*> final_hists;
    final_hists.push_back(std::move(new_data_hist));
    final_hists.push_back(std::move(new_mc_hist));

    return final_hists;    
}

//=============================================================================================================================

void plot_all_planes(string X, string Y, bool profile, int& argc, char** argv) {
    string filename1 = "histogramsDATA.root";
    string filename2 = "histogramsMC.root";

    TFile* file1 = TFile::Open(filename1.c_str(), "READ");
    TFile* file2 = TFile::Open(filename2.c_str(), "READ");

    auto data_histograms = loop_through_root_file(file1, X, Y);
    auto mc_histograms = loop_through_root_file(file2, X, Y);

    unique_ptr<TH2D> plane0data = nullptr;
    unique_ptr<TH2D> plane1data = nullptr;
    unique_ptr<TH2D> plane2data = nullptr;

    combineHistogramsWithName(data_histograms, "0", plane0data);
    combineHistogramsWithName(data_histograms, "1", plane1data);
    combineHistogramsWithName(data_histograms, "2", plane2data);

    unique_ptr<TH2D> plane0mc = nullptr;
    unique_ptr<TH2D> plane1mc = nullptr;
    unique_ptr<TH2D> plane2mc = nullptr;

    combineHistogramsWithName(mc_histograms, "0", plane0mc);
    combineHistogramsWithName(mc_histograms, "1", plane1mc);
    combineHistogramsWithName(mc_histograms, "2", plane2mc);

    // I find each plane has a sweet spot for the minimum data found thorugh trial and error
    //
    auto plane0hists = adjust_bins(plane0data,plane0mc, 500.);
    auto plane1hists = adjust_bins(plane1data,plane1mc, 400.);
    auto plane2hists = adjust_bins(plane2data,plane2mc, 300.);

    auto profile0d = plane0hists[0]->ProfileX(("pfx01"));
    auto profile1d = plane1hists[0]->ProfileX(("pfx11"));
    auto profile2d = plane2hists[0]->ProfileX(("pfx21"));

    auto profile0m = plane0hists[1]->ProfileX(("pfx0"));
    auto profile1m = plane1hists[1]->ProfileX(("pfx1"));
    auto profile2m = plane2hists[1]->ProfileX(("pfx2"));

    auto profile0 = ratio_profile(profile0d,profile0m);
    profile0->SetName("plane 0");
    auto profile1 = ratio_profile(profile1d,profile1m);
    profile1->SetName("plane 1");
    auto profile2 = ratio_profile(profile2d,profile2m);
    profile2->SetName("plane 2");



    //set up titles for plots and remove stat box
    string title0 = X + " vs. " + Y + " Plane 0 Ratio (Data/MC)";
    string title1 = X + " vs. " + Y + " Plane 1 Ratio (Data/MC)";
    string title2 = X + " vs. " + Y + " Plane 2 Ratio (Data/MC)";
    string xAxis;
    if (X == "Theta" || X== "Phi"){
        xAxis = X + " (radians)";
    }
    else{
        xAxis = X + " (cm)";
    }
    string yAxis = Y;
    profile0->SetTitle(title0.c_str());
    profile1->SetTitle(title1.c_str());
    profile2->SetTitle(title2.c_str());
    profile0->SetXTitle(xAxis.c_str());
    profile1->SetXTitle(xAxis.c_str());
    profile2->SetXTitle(xAxis.c_str());
    profile0->SetYTitle(yAxis.c_str());
    profile1->SetYTitle(yAxis.c_str());
    profile2->SetYTitle(yAxis.c_str());
    profile0->SetStats(false);
    profile1->SetStats(false);
    profile2->SetStats(false);

    plane0hists[0]->Divide(plane0hists[1]);
    plane1hists[0]->Divide(plane1hists[1]);
    plane2hists[0]->Divide(plane2hists[1]);

    plane0hists[0]->SetTitle(title0.c_str());
    plane1hists[0]->SetTitle(title1.c_str());
    plane2hists[0]->SetTitle(title2.c_str());
    plane0hists[0]->SetXTitle(xAxis.c_str());
    plane1hists[0]->SetXTitle(xAxis.c_str());
    plane2hists[0]->SetXTitle(xAxis.c_str());
    plane0hists[0]->SetYTitle(yAxis.c_str());
    plane1hists[0]->SetYTitle(yAxis.c_str());
    plane2hists[0]->SetYTitle(yAxis.c_str());
    plane0hists[0]->SetStats(false);
    plane1hists[0]->SetStats(false);
    plane2hists[0]->SetStats(false);

    //set ranges for histograms then profiles to give same color bar or just set max
    double global_max_hist = std::min(std::min(plane0hists[0]->GetMaximum(),plane1hists[0]->GetMaximum()),plane2hists[0]->GetMaximum());
    double global_min_hist = std::min(std::min(plane0hists[0]->GetMinimum(),plane1hists[0]->GetMinimum()),plane2hists[0]->GetMinimum());
    plane0hists[0]->SetMaximum(global_max_hist);
    plane1hists[0]->SetMaximum(global_max_hist);
    plane2hists[0]->SetMaximum(global_max_hist);
    plane0hists[0]->SetMinimum(global_min_hist);
    plane1hists[0]->SetMinimum(global_min_hist);
    plane2hists[0]->SetMinimum(global_min_hist);

    /*
    double global_max_prof = std::max(std::max(profile0->GetMaximum(),profile1->GetMaximum()),profile2->GetMaximum());
    double global_min_prof = std::min(std::min(profile0->GetMinimum(),profile1->GetMinimum()),profile2->GetMinimum());
    profile0->SetMaximum(global_max_prof);
    profile1->SetMaximum(global_max_prof);
    profile2->SetMaximum(global_max_prof);
    profile0->SetMinimum(global_min_prof);
    profile1->SetMinimum(global_min_prof);
    profile2->SetMinimum(global_min_prof);
    */


    TApplication app("app", &argc, argv);

    auto canvas = new TCanvas("canvas", "Histograms", 800, 600);
    if (profile){
        canvas->Divide(1,3);
        canvas->cd(1); profile0->Draw("COLZ");
        canvas->cd(2); profile1->Draw("COLZ");
        canvas->cd(3); profile2->Draw("COLZ");
        canvas->Update();
    }
    else{
        canvas->Divide(1,3);
        canvas->cd(1); plane0hists[0]->Draw("COLZ");
        canvas->cd(2); plane1hists[0]->Draw("COLZ");
        canvas->cd(3); plane2hists[0]->Draw("COLZ");
        canvas->Update();
        }

    app.Run();

    delete canvas;

    file1->Close();
    file2->Close();

    delete file2;
    delete file1;
    }

//=============================================================================================================================

// Function to redistribute counts from an old 2D histogram to a new binning scheme using TArrayD
TH3D* RedistributeCounts3D(TH3D* old_hist, TArrayD* x_edges_new, TArrayD* y_edges_new, TArrayD* z_edges_new) {
    // Create a new histogram with the new binning scheme
    TH3D* new_hist =  new TH3D("new_hist", "Redistributed Histogram",
                                           x_edges_new->GetSize() - 1, x_edges_new->GetArray(),
                                           y_edges_new->GetSize() - 1, y_edges_new->GetArray(),
                                           z_edges_new->GetSize() - 1, z_edges_new->GetArray());

    // Iterate through all old bins
    for (int i = 1; i <= old_hist->GetNbinsX(); ++i) {
        for (int j = 1; j <= old_hist->GetNbinsY(); ++j) {
            for (int k = 1; k <= old_hist->GetNbinsZ(); ++k){
                double x_center = old_hist->GetXaxis()->GetBinCenter(i);
                double y_center = old_hist->GetYaxis()->GetBinCenter(j);
                double z_center = old_hist->GetZaxis()->GetBinCenter(k);
                double content = old_hist->GetBinContent(i, j, k);

                // Assuming the old bins map directly to new bins or fall entirely within new bins
                // Find the new bin indices for the given bin center
                int new_i = new_hist->GetXaxis()->FindBin(x_center);
                int new_j = new_hist->GetYaxis()->FindBin(y_center);
                int new_k = new_hist->GetZaxis()->FindBin(z_center);

                // Check if the bin center is within the bounds of the new histogram
                if (new_i > 0 && new_i <= new_hist->GetNbinsX() && new_j > 0 && new_j <= new_hist->GetNbinsY() && new_k > 0 && new_k <= new_hist->GetNbinsZ()) {
                    // Add the old bin's content to the corresponding new bin
                    double existingContent = new_hist->GetBinContent(new_i, new_j,new_k);
                    new_hist->SetBinContent(new_i, new_j,new_k, existingContent + content);
                }
            }
        }
    }
    return new_hist;
}

//=============================================================================================================================

//function for taking a histogram and adjusting the bins so that there are no entries with zero
vector<TArrayD*> remove_zeros3D( TH3D* old_hist1, TH3D* old_hist2) {
    // Access bin arrays
    int nxbins = old_hist1->GetXaxis()->GetNbins();
    int nybins = old_hist1->GetYaxis()->GetNbins();
    int nzbins = old_hist1->GetZaxis()->GetNbins();

    //make vector of bins because they are easy to work with
    vector<double> xbins, ybins, zbins;

    // Fill the vectors with the bin edges
    for (int i = 0; i <= nxbins - 1; i++) {
        xbins.push_back(old_hist1->GetXaxis()->GetBinLowEdge(i+1));
    }

    for (int j = 0; j <= nybins - 1; j++) {
        ybins.push_back(old_hist1->GetYaxis()->GetBinLowEdge(j+1));
    }

    for (int k = 0; k <= nzbins - 1; k++) {
        zbins.push_back(old_hist1->GetZaxis()->GetBinLowEdge(k+1));
    }

    // Clone the original histogram for manipulation
    TH3D* new_hist1 = new TH3D(*old_hist1);
    TH3D* new_hist2 = new TH3D(*old_hist2);

    bool zeros1 = hasZeroBinEntries3D(old_hist1);
    bool zeros2 = hasZeroBinEntries3D(old_hist2);

    // Create vector that will contain everytime there is a bin with zero in it, then the x and y bins with the most amount of zeros will be removed 
    vector<int> xzeros, yzeros, zzeros;

    //add count to stop infinte loop
    int count = 0;

    //loops thorugh histogram and removes most the x and y bins with the most zeros until no zeros remain (or until nxbins is reached because this loop can go on forever and that seems like a good cutoff)
    while (zeros1 || zeros2){
        count += 1;
        //Make a vector with an ent
        for (int i = 0 ; i <= xbins.size(); i++){ 
            for (int j = 0; j <= ybins.size() ; j++){
                for (int k = 0; k <= zbins.size() ; k++){
                    Double_t content1 = new_hist1->GetBinContent(i + 1, j + 1, k + 1);
                    Double_t content2 = new_hist2->GetBinContent(i + 1, j + 1, k + 1);
                    Double_t content = content1 + content2;
                    if (content == 0){
                        xzeros.push_back(i);
                        yzeros.push_back(j);
                        zzeros.push_back(k);
                    }
                }
            }
        }
        int i = findMode(xzeros);
        int j = findMode(yzeros);
        int k = findMode(zzeros);
        //remove bins with the most zeros if they are not too small (larger here since these have more points)
        if (xbins.size() > 10){
            xbins.erase(xbins.begin() + i - 1);
        }
        if (ybins.size() > 10){
            ybins.erase(ybins.begin()+ j - 1);
        }
        if (zbins.size() > 10){
            zbins.erase(zbins.begin()+ k - 1);
        }

         // Creating a new histogram to renter the loop with the bins removed
        //make a new array 
        TArrayD* xbin_array = new TArrayD(xbins.size(), xbins.data());
        TArrayD* ybin_array = new TArrayD(ybins.size(), ybins.data());
        TArrayD* zbin_array = new TArrayD(zbins.size(), zbins.data());

        new_hist1 = RedistributeCounts3D(new_hist1, xbin_array, ybin_array, zbin_array);
        new_hist2 = RedistributeCounts3D(new_hist2, xbin_array, ybin_array, zbin_array);

        delete xbin_array;
        delete ybin_array;
        delete zbin_array;

        xzeros.clear();
        yzeros.clear();
        zzeros.clear();

        zeros1 = hasZeroBinEntries3D(new_hist1);
        zeros2 = hasZeroBinEntries3D(new_hist2);

        //stop code from running forever if zeros are not removed
        if (count > nxbins){
            cout<< "Can't Remove all Zeros" << endl;
            break;
        }
    }
    //return bins 
    TArrayD* xbin_array = new TArrayD(xbins.size(), xbins.data());
    TArrayD* ybin_array = new TArrayD(ybins.size(), ybins.data());
    TArrayD* zbin_array = new TArrayD(zbins.size(), zbins.data());

    vector<TArrayD*> bins = {xbin_array,ybin_array, zbin_array};
    return bins;
}

//=============================================================================================================================

//adjust the bins for the 3D histograms same logic as 2D version but with 1 more of everyting essentially
void adjust_bins3D( unique_ptr<TH3D>& data_hist, unique_ptr<TH3D>& mc_hist, const double& min_data){
     //----------------------------------------------------------------
    // create bins for histograms that don't have array of them
    int nBinsX = data_hist->GetXaxis()->GetNbins();
    double xMin = data_hist->GetXaxis()->GetXmin();
    double xMax = data_hist->GetXaxis()->GetXmax();

    int nBinsY = data_hist->GetYaxis()->GetNbins();
    double yMin = data_hist->GetYaxis()->GetXmin();
    double yMax = data_hist->GetYaxis()->GetXmax();

    int nBinsZ = data_hist->GetZaxis()->GetNbins();
    double zMin = data_hist->GetZaxis()->GetXmin();
    double zMax = data_hist->GetZaxis()->GetXmax();
  
    TArrayD* xbins = CreateUniformBinArray(nBinsX, xMin, xMax);
    TArrayD* ybins = CreateUniformBinArray(nBinsY, yMin, yMax);
    TArrayD* zbins = CreateUniformBinArray(nBinsZ, zMin, zMax);
     //----------------------------------------------------------------
    // adjusts the x and y bins for the minimum amount of data needed
    TH1D* data_x_proj = data_hist->ProjectionX();
    TH1D* data_y_proj = data_hist->ProjectionY();
    TH1D* data_z_proj = data_hist->ProjectionZ();
    TH1D* mc_x_proj = mc_hist->ProjectionX();
    TH1D* mc_y_proj = mc_hist->ProjectionY();
    TH1D* mc_z_proj = mc_hist->ProjectionZ();

    auto new_xbins = create_new_bins(xbins,data_x_proj,mc_x_proj,min_data);
    auto new_ybins = create_new_bins(ybins,data_y_proj,mc_y_proj,min_data);
    auto new_zbins = create_new_bins(zbins,data_z_proj,mc_z_proj,min_data);

    TH3D* data_hist_var = data_hist.get();
    TH3D* mc_hist_var = mc_hist.get();

    // Use the redistributed counts to create new histograms
    auto new_data_hist = RedistributeCounts3D(data_hist_var, new_xbins, new_ybins,new_zbins);
    auto new_mc_hist = RedistributeCounts3D(mc_hist_var, new_xbins, new_ybins,new_zbins);
    // ----------------------------------------------------------------

    //remove the zeros from the histogram
    auto newest_bins = remove_zeros3D(new_data_hist,new_mc_hist);
    
    new_data_hist = RedistributeCounts3D(new_data_hist, newest_bins[0], newest_bins[1], newest_bins[2]);
    new_mc_hist = RedistributeCounts3D(new_mc_hist, newest_bins[0], newest_bins[1], newest_bins[2]);

    // ----------------------------------------------------------------
    // change inputs to new histogram
    data_hist.reset(new_data_hist);
    mc_hist.reset(new_mc_hist);
}

//=============================================================================================================================

void plot_all_planes_3D(string X, string Y, string Z, bool profile,int& argc, char** argv){

    string filename1 = "histogramsDATA.root";
    string filename2 = "histogramsMC.root";

    TFile* file1 = TFile::Open(filename1.c_str(), "READ");
    TFile* file2 = TFile::Open(filename2.c_str(), "READ");

     //-----------------------------------------------------
    // set up 2D histograms
    auto data_XZ = loop_through_root_file(file1, X, Z);
    auto mc_XZ = loop_through_root_file(file2, X, Z);
    auto data_YZ = loop_through_root_file(file1, Y, Z);
    auto mc_YZ = loop_through_root_file(file2, Y, Z);


    unique_ptr<TH2D> plane0data_xz = nullptr;
    unique_ptr<TH2D> plane1data_xz = nullptr;
    unique_ptr<TH2D> plane2data_xz = nullptr;
    combineHistogramsWithName(data_XZ, "0", plane0data_xz);
    combineHistogramsWithName(data_XZ, "1", plane1data_xz);
    combineHistogramsWithName(data_XZ, "2", plane2data_xz);

    unique_ptr<TH2D> plane0mc_xz = nullptr;
    unique_ptr<TH2D> plane1mc_xz = nullptr;
    unique_ptr<TH2D> plane2mc_xz = nullptr;
    combineHistogramsWithName(mc_XZ, "0", plane0mc_xz);
    combineHistogramsWithName(mc_XZ, "1", plane1mc_xz);
    combineHistogramsWithName(mc_XZ, "2", plane2mc_xz);

    unique_ptr<TH2D> plane0data_yz = nullptr;
    unique_ptr<TH2D> plane1data_yz = nullptr;
    unique_ptr<TH2D> plane2data_yz = nullptr;
    combineHistogramsWithName(data_YZ, "0", plane0data_yz);
    combineHistogramsWithName(data_YZ, "1", plane1data_yz);
    combineHistogramsWithName(data_YZ, "2", plane2data_yz);

    unique_ptr<TH2D> plane0mc_yz = nullptr;
    unique_ptr<TH2D> plane1mc_yz = nullptr;
    unique_ptr<TH2D> plane2mc_yz = nullptr;
    combineHistogramsWithName(mc_YZ , "0", plane0mc_yz);
    combineHistogramsWithName(mc_YZ , "1", plane1mc_yz);
    combineHistogramsWithName(mc_YZ , "2", plane2mc_yz);

    //---------------------------------------------------------------
    // Try not adjusting bins can change
    // set up 1D hists
    unique_ptr<TH1D> plane0data_x(plane0data_xz->ProjectionX("plane0data_x"));
    unique_ptr<TH1D> plane0data_z1(plane0data_xz->ProjectionY("plane0data_z1"));
    unique_ptr<TH1D> plane1data_x(plane1data_xz->ProjectionX("plane1data_x"));
    unique_ptr<TH1D> plane1data_z1(plane1data_xz->ProjectionY("plane1data_z1"));
    unique_ptr<TH1D> plane2data_x(plane2data_xz->ProjectionX("plane2data_x"));
    unique_ptr<TH1D> plane2data_z1(plane2data_xz->ProjectionY("plane2data_z1"));
    unique_ptr<TH1D> plane0data_y(plane0data_yz->ProjectionX("plane0data_y"));
    unique_ptr<TH1D> plane0data_z2(plane0data_yz->ProjectionY("plane0data_z2"));
    unique_ptr<TH1D> plane1data_y(plane1data_yz->ProjectionX("plane1data_y"));
    unique_ptr<TH1D> plane1data_z2(plane1data_yz->ProjectionY("plane1data_z2"));
    unique_ptr<TH1D> plane2data_y(plane2data_yz->ProjectionX("plane2data_y"));
    unique_ptr<TH1D> plane2data_z2(plane2data_yz->ProjectionY("plane2data_z2"));
    unique_ptr<TH1D> plane0mc_x(plane0mc_xz->ProjectionX("plane0mc_x"));
    unique_ptr<TH1D> plane0mc_z1(plane0mc_xz->ProjectionY("plane0mc_z1"));
    unique_ptr<TH1D> plane1mc_x(plane1mc_xz->ProjectionX("plane1mc_x"));
    unique_ptr<TH1D> plane1mc_z1(plane1mc_xz->ProjectionY("plane1mc_z1"));
    unique_ptr<TH1D> plane2mc_x(plane2mc_xz->ProjectionX("plane2mc_x"));
    unique_ptr<TH1D> plane2mc_z1(plane2mc_xz->ProjectionY("plane2mc_z1"));
    unique_ptr<TH1D> plane0mc_y(plane0mc_yz->ProjectionX("plane0mc_y"));
    unique_ptr<TH1D> plane0mc_z2(plane0mc_yz->ProjectionY("plane0mc_z2"));
    unique_ptr<TH1D> plane1mc_y(plane1mc_yz->ProjectionX("plane1mc_y"));
    unique_ptr<TH1D> plane1mc_z2(plane1mc_yz->ProjectionY("plane1mc_z2"));
    unique_ptr<TH1D> plane2mc_y(plane2mc_yz->ProjectionX("plane2mc_y"));
    unique_ptr<TH1D> plane2mc_z2(plane2mc_yz->ProjectionY("plane2mc_z2"));
    //----------------------------------------------------------------
    //make 3D hists
    unique_ptr<TH3D> plane0_3D_data = make_unique<TH3D>("Plane 0 X Y X Data", "Plane 0 Data", 
                                        plane0data_x->GetXaxis()->GetNbins(), plane0data_x->GetXaxis()->GetXmin(), plane0data_x->GetXaxis()->GetXmax(), 
                                        plane0data_y->GetXaxis()->GetNbins(), plane0data_y->GetXaxis()->GetXmin(), plane0data_y->GetXaxis()->GetXmax(), 
                                        plane0data_z1->GetXaxis()->GetNbins(), plane0data_z1->GetXaxis()->GetXmin(), plane0data_z1->GetXaxis()->GetXmax() );
    for (int i = 1; i <= plane0data_x->GetXaxis()->GetNbins(); i++){
        for (int j = 1; j <= plane0data_y->GetXaxis()->GetNbins(); j++){
            for (int k = 1; k <= plane0data_z1->GetXaxis()->GetNbins(); k++){
                plane0_3D_data->Fill(
                plane0data_x->GetXaxis()->GetBinCenter(i),
                plane0data_y->GetXaxis()->GetBinCenter(j),
                plane0data_z1->GetXaxis()->GetBinCenter(k),
                plane0data_x->GetBinContent(i) + plane0data_y->GetBinContent(j) + plane0data_z1->GetBinContent(k) + plane0data_z2->GetBinContent(k)
                );
            }
        }
    };
    unique_ptr<TH3D> plane1_3D_data = make_unique<TH3D>("Plane 1 X Y X Data", "Plane 1 Data", 
                                        plane1data_x->GetXaxis()->GetNbins(), plane1data_x->GetXaxis()->GetXmin(), plane1data_x->GetXaxis()->GetXmax(), 
                                        plane1data_y->GetXaxis()->GetNbins(), plane1data_y->GetXaxis()->GetXmin(), plane1data_y->GetXaxis()->GetXmax(), 
                                        plane1data_z1->GetXaxis()->GetNbins(), plane1data_z1->GetXaxis()->GetXmin(), plane1data_z1->GetXaxis()->GetXmax() );
    for (int i = 1; i <= plane1data_x->GetXaxis()->GetNbins(); i++){
        for (int j = 1; j <= plane1data_y->GetXaxis()->GetNbins(); j++){
            for (int k = 1; k <= plane1data_z1->GetXaxis()->GetNbins(); k++){
                plane1_3D_data->Fill(
                plane1data_x->GetXaxis()->GetBinCenter(i),
                plane1data_y->GetXaxis()->GetBinCenter(j),
                plane1data_z1->GetXaxis()->GetBinCenter(k),
                plane1data_x->GetBinContent(i) + plane1data_y->GetBinContent(j) + plane1data_z1->GetBinContent(k) + plane1data_z2->GetBinContent(k)
                );
            }
        }
    };
    unique_ptr<TH3D> plane2_3D_data = make_unique<TH3D>("Plane 2 X Y X Data", "Plane 2 Data", 
                                        plane2data_x->GetXaxis()->GetNbins(), plane2data_x->GetXaxis()->GetXmin(), plane2data_x->GetXaxis()->GetXmax(), 
                                        plane2data_y->GetXaxis()->GetNbins(), plane2data_y->GetXaxis()->GetXmin(), plane2data_y->GetXaxis()->GetXmax(), 
                                        plane2data_z1->GetXaxis()->GetNbins(), plane2data_z1->GetXaxis()->GetXmin(), plane2data_z1->GetXaxis()->GetXmax() );
    for (int i = 1; i <= plane2data_x->GetXaxis()->GetNbins(); i++){
        for (int j = 1; j <= plane2data_y->GetXaxis()->GetNbins(); j++){
            for (int k = 1; k <= plane2data_z1->GetXaxis()->GetNbins(); k++){
                plane2_3D_data->Fill(
                plane2data_x->GetXaxis()->GetBinCenter(i),
                plane2data_y->GetXaxis()->GetBinCenter(j),
                plane2data_z1->GetXaxis()->GetBinCenter(k),
                plane2data_x->GetBinContent(i) + plane2data_y->GetBinContent(j) + plane2data_z1->GetBinContent(k) + plane2data_z2->GetBinContent(k)
                );
            }
        }
    };
    unique_ptr<TH3D> plane0_3D_mc = make_unique<TH3D>("Plane 0 X Y X MC", "Plane 0 MC", 
                                        plane0mc_x->GetXaxis()->GetNbins(), plane0mc_x->GetXaxis()->GetXmin(), plane0mc_x->GetXaxis()->GetXmax(), 
                                        plane0mc_y->GetXaxis()->GetNbins(), plane0mc_y->GetXaxis()->GetXmin(), plane0mc_y->GetXaxis()->GetXmax(), 
                                        plane0mc_z1->GetXaxis()->GetNbins(), plane0mc_z1->GetXaxis()->GetXmin(), plane0mc_z1->GetXaxis()->GetXmax() );
    for (int i = 1; i <= plane0mc_x->GetXaxis()->GetNbins(); i++){
        for (int j = 1; j <= plane0mc_y->GetXaxis()->GetNbins(); j++){
            for (int k = 1; k <= plane0mc_z1->GetXaxis()->GetNbins(); k++){
                plane0_3D_mc->Fill(
                plane0mc_x->GetXaxis()->GetBinCenter(i),
                plane0mc_y->GetXaxis()->GetBinCenter(j),
                plane0mc_z1->GetXaxis()->GetBinCenter(k),
                plane0mc_x->GetBinContent(i) + plane0mc_y->GetBinContent(j) + plane0mc_z1->GetBinContent(k) + plane0mc_z2->GetBinContent(k)
                );
            }
        }
    };
    unique_ptr<TH3D> plane1_3D_mc = make_unique<TH3D>("Plane 1 X Y X MC", "Plane 1 MC", 
                                        plane1mc_x->GetXaxis()->GetNbins(), plane1mc_x->GetXaxis()->GetXmin(), plane1mc_x->GetXaxis()->GetXmax(), 
                                        plane1mc_y->GetXaxis()->GetNbins(), plane1mc_y->GetXaxis()->GetXmin(), plane1mc_y->GetXaxis()->GetXmax(), 
                                        plane1mc_z1->GetXaxis()->GetNbins(), plane1mc_z1->GetXaxis()->GetXmin(), plane1mc_z1->GetXaxis()->GetXmax() );
    for (int i = 1; i <= plane1mc_x->GetXaxis()->GetNbins(); i++){
        for (int j = 1; j <= plane1mc_y->GetXaxis()->GetNbins(); j++){
            for (int k = 1; k <= plane1mc_z1->GetXaxis()->GetNbins(); k++){
                plane1_3D_mc->Fill(
                plane1mc_x->GetXaxis()->GetBinCenter(i),
                plane1mc_y->GetXaxis()->GetBinCenter(j),
                plane1mc_z1->GetXaxis()->GetBinCenter(k),
                plane1mc_x->GetBinContent(i) + plane1mc_y->GetBinContent(j) + plane1mc_z1->GetBinContent(k) + plane1mc_z2->GetBinContent(k)
                );
            }
        }
    };
    unique_ptr<TH3D> plane2_3D_mc = make_unique<TH3D>("Plane 2 X Y X MC", "Plane 2 MC", 
                                        plane2mc_x->GetXaxis()->GetNbins(), plane2mc_x->GetXaxis()->GetXmin(), plane2mc_x->GetXaxis()->GetXmax(), 
                                        plane2mc_y->GetXaxis()->GetNbins(), plane2mc_y->GetXaxis()->GetXmin(), plane2mc_y->GetXaxis()->GetXmax(), 
                                        plane2mc_z1->GetXaxis()->GetNbins(), plane2mc_z1->GetXaxis()->GetXmin(), plane2mc_z1->GetXaxis()->GetXmax() );
    for (int i = 1; i <= plane2mc_x->GetXaxis()->GetNbins(); i++){
        for (int j = 1; j <= plane2mc_y->GetXaxis()->GetNbins(); j++){
            for (int k = 1; k <= plane2mc_z1->GetXaxis()->GetNbins(); k++){
                plane2_3D_mc->Fill(
                plane2mc_x->GetXaxis()->GetBinCenter(i),
                plane2mc_y->GetXaxis()->GetBinCenter(j),
                plane2mc_z1->GetXaxis()->GetBinCenter(k),
                plane2mc_x->GetBinContent(i) + plane2mc_y->GetBinContent(j) + plane2mc_z1->GetBinContent(k) + plane2mc_z2->GetBinContent(k)
                );
            }
        }
    };
    adjust_bins3D(plane0_3D_data,plane0_3D_mc, 1000.);
    adjust_bins3D(plane1_3D_data,plane1_3D_mc, 800.);
    adjust_bins3D(plane2_3D_data,plane2_3D_mc, 600.);
     //-----------------------------------------------------------------------
    //plot either hist 3d ratios or profile ratios
    string title0 = X + " vs. " + Y + " vs. " + Z  + " Plane 0 Ratio (Data/MC)";
    string title1 = X + " vs. " + Y + " vs. " + Z  + " Plane 1 Ratio (Data/MC)";
    string title2 = X + " vs. " + Y + " vs. " + Z  + " Plane 2 Ratio (Data/MC)";
    string zAxis = Z;
    string xAxis;
    string yAxis;
    if (X == "Theta" || X== "Phi"){
        xAxis = X + " (radians)";
    }
    else{
        xAxis = X + " (cm)";
    }
    if (Y == "Theta" || Y== "Phi"){
        yAxis = Y + " (radians)";
    }
    else{
        yAxis = Y + " (cm)";
    }

    if (profile){
        unique_ptr<TProfile2D> profile0(plane0_3D_data->Project3DProfile("xy"));
        profile0->SetName("Plane 0");
        unique_ptr<TProfile2D> profile00(plane0_3D_mc->Project3DProfile("xy"));
        profile00->SetName("Plane 0 mc");
        unique_ptr<TProfile2D> profile1(plane1_3D_data->Project3DProfile("xy"));
        profile1->SetName("Plane 1");
        unique_ptr<TProfile2D> profile11(plane1_3D_mc->Project3DProfile("xy"));
        profile11->SetName("Plane 1 mc");
        unique_ptr<TProfile2D> profile2(plane2_3D_data->Project3DProfile("xy"));
        profile2->SetName("Plane 2");
        unique_ptr<TProfile2D> profile22(plane2_3D_mc->Project3DProfile("xy"));
        profile22->SetName("Plane 2 mc");

        profile0->Divide(profile00.get());
        profile1->Divide(profile11.get());
        profile2->Divide(profile22.get());

        profile0->SetTitle(title0.c_str());
        profile1->SetTitle(title1.c_str());
        profile2->SetTitle(title2.c_str());
        profile0->SetXTitle(xAxis.c_str());
        profile1->SetXTitle(xAxis.c_str());
        profile2->SetXTitle(xAxis.c_str());
        profile0->SetYTitle(yAxis.c_str());
        profile1->SetYTitle(yAxis.c_str());
        profile2->SetYTitle(yAxis.c_str());
        profile0->SetStats(false);
        profile1->SetStats(false);
        profile2->SetStats(false);

        //set ranges for profiles to give same color bar
        double global_max_prof = std::max(std::max(profile0->GetMaximum(),profile1->GetMaximum()),profile2->GetMaximum());
        double global_min_prof = std::min(std::min(profile0->GetMinimum(),profile1->GetMinimum()),profile2->GetMinimum());
        profile0->SetMaximum(global_max_prof);
        profile1->SetMaximum(global_max_prof);
        profile2->SetMaximum(global_max_prof);
        profile0->SetMinimum(global_min_prof);
        profile1->SetMinimum(global_min_prof);
        profile2->SetMinimum(global_min_prof);

        TApplication app("app", &argc, argv);

        auto canvas = new TCanvas("canvas", "Histograms", 800, 600);

        canvas->Divide(1,3);
        canvas->cd(1); profile0->Draw("COLZ");
        canvas->cd(2); profile1->Draw("COLZ");
        canvas->cd(3); profile2->Draw("COLZ");
        canvas->Update();

        app.Run();

        delete canvas;
    }
    else{
        plane0_3D_data->Divide(plane0_3D_mc.get());
        plane1_3D_data->Divide(plane1_3D_mc.get());
        plane2_3D_data->Divide(plane2_3D_mc.get());

        plane0_3D_data->SetTitle(title0.c_str());
        plane1_3D_data->SetTitle(title1.c_str());
        plane2_3D_data->SetTitle(title2.c_str());
        plane0_3D_data->SetXTitle(xAxis.c_str());
        plane1_3D_data->SetXTitle(xAxis.c_str());
        plane2_3D_data->SetXTitle(xAxis.c_str());
        plane0_3D_data->SetYTitle(yAxis.c_str());
        plane1_3D_data->SetYTitle(yAxis.c_str());
        plane2_3D_data->SetYTitle(yAxis.c_str());
        plane0_3D_data->SetZTitle(zAxis.c_str());
        plane1_3D_data->SetZTitle(zAxis.c_str());
        plane2_3D_data->SetZTitle(zAxis.c_str());
        plane0_3D_data->SetStats(false);
        plane1_3D_data->SetStats(false);
        plane2_3D_data->SetStats(false);

        TApplication app("app", &argc, argv);

        auto canvas = new TCanvas("canvas", "Histograms", 800, 600);

        canvas->Divide(1,3);
        canvas->cd(1); plane0_3D_data->Draw("COLZ");
        canvas->cd(2); plane1_3D_data->Draw("COLZ");
        canvas->cd(3); plane2_3D_data->Draw("COLZ");
        canvas->Update();

        app.Run();

        delete canvas;
    }

    file1->Close();
    file2->Close();

    delete file2;
    delete file1;
}