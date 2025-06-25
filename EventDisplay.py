import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from matplotlib.patches import Circle
from matplotlib.collections import PatchCollection
import matplotlib.colors as colors
import copy
from contextlib import nullcontext
import os

class EventDisplay:
  
    def load_mPMT_positions(self,fileName):   
        
        # Get the directory where this script (module) is located
        module_dir = os.path.dirname(os.path.abspath(__file__))
        
        # Construct the full path to the file
        file_path = os.path.join(module_dir, fileName)
        
        #load the mpmt_positions file and create a series of variables used throughout the 
        self.mPMT_2D_projection = np.loadtxt(file_path, delimiter=',', skiprows=1, dtype = int)
        self.nmPMTs = len(self.mPMT_2D_projection[:,0])

        self.nChannels = 19*self.nmPMTs
        self.mask_list = None
    
    def load_wcsim_tubeno_mapping(self,fileName):
        # Get the directory where this script (module) is located
        module_dir = os.path.dirname(os.path.abspath(__file__))
        
        # Construct the full path to the file
        file_path = os.path.join(module_dir, fileName)
        
        data = np.loadtxt(file_path, skiprows=5, usecols=(0, 1, 2),dtype=int)
        self.wcsim_tube_mapping = {row[0]: (row[1], row[2]) for row in data}
        
        # #load the mpmt_positions file and create a series of variables used throughout the 
        # self.mPMT_2D_projection = np.loadtxt(file_path, delimiter=',', skiprows=1, dtype = int)
        # self.nmPMTs = len(self.mPMT_2D_projection[:,0])

        # self.nChannels = 19*self.nmPMTs
        # self.mask_list = None
    
    def map_wcsim_tubeno_to_slot_pmt_id(self, tube_nos):
        
        if(self.wcsim_tube_mapping == None):
            raise Exception("wcsim_tube_mapping not loaded in and calling mapping")
        
        # mPMT_id, PMT_pos = map(list,zip(*[self.wcsim_tube_mapping[tube_no] for tube_no in tube_nos]))
        mPMT_id, PMT_pos = np.array([self.wcsim_tube_mapping[tube_no] for tube_no in tube_nos]).T
        print(mPMT_id)
        #in the geo file PMT pos is defined starting at 1 to -1 to start from 0
        PMT_pos = PMT_pos-1
        return mPMT_id, PMT_pos

            
    def mask_mPMTs(self,mask_list):
        self.mask_list = mask_list

    def process_data(self,pmts, pmt_data,sum_data=True,average_data=False):
        #takes as input a list of pmts and pmt_data
        #returns a vector of length n_mpmts*19 with the charge in each column
        #if the sum_data is true then if we have multiple of the same channel in pmts then the entries are summed for that channel
        
        outData = np.zeros(self.nChannels)
        hitCount = np.zeros(self.nChannels)

        for iPMT, channel_data in zip(pmts,pmt_data):
            #print(channel_data)
            if(channel_data > 0):
                #print(channel_data, hitCount[iPMT], iPMT)
                hitCount[iPMT] += 1

            if(sum_data): 
                if (channel_data > 0):
                    outData[iPMT]+=channel_data
            else:
                #chose the smallest value of data for that channel e.g time  
                if(outData[iPMT]>0):
                    outData[iPMT]=min(outData[iPMT],channel_data)
                else:
                    outData[iPMT]=channel_data

        if average_data:
            for index in range(len(hitCount)):
                if hitCount[index] > 1:
                    outData[index] = outData[index]/hitCount[index]
                else:
                    outData[index] = 0

#            outData /= len(pmt_data)            
     
        return outData

    def channel_position_offset(self, channel, rotationAngle):
        """
        Calculate offset in plotting coordinates for channel of PMT within mPMT.

        Parameters
        ----------
        channel: array_like of float
            array of channel IDs or PMT IDs
        use_new_convention: bool
            use newer convention for the channel mapping of PMTs in the mPMT (starts with central PMT, then middle ring then
            outer ring) as opposed to old convention (starts with outer ring, then middle ring, then central PMT).

        Returns
        -------
        np.ndarray:
            Array of (x, y) coordinate offsets
        """
        channel = channel % 19
        theta = (channel > 6)*2*np.pi*(19-channel)/12 + ((channel > 0) & (channel <= 6))*2*np.pi*(7-channel)/6 - rotationAngle*np.pi/180
        radius = 0.2*(channel > 0) + 0.2*(channel > 6)

        position = np.column_stack((radius*np.sin(theta), radius*np.cos(theta)))
        return position

    def coordinates_eachChannel(self):
        #returns an array of size coordinates[nChannels,2] with the 2d position on the event display for each channel
        coordinates = np.zeros((self.nChannels,2))        
        #vectors of length self.nChannels, for each channel get the channel no. mPMT number, rotation angle of the mPMT 
        channels = np.array([19*mPMT+j for mPMT in self.mPMT_2D_projection[:,0] for j in range(19)])
        mPMTs = np.array([mPMT for mPMT in self.mPMT_2D_projection[:,0] for j in range(19)])
        rotationAngle = np.array([mPMT for mPMT in self.mPMT_2D_projection[:,3] for j in range(19)])
        
        #sets the coords of each channel at the centre of the mPMT
        coordinates[:,0] = [self.mPMT_2D_projection[mPMT,1] for mPMT in mPMTs]
        coordinates[:,1] = [self.mPMT_2D_projection[mPMT,2] for mPMT in mPMTs]

        coordinates += self.channel_position_offset(channels,rotationAngle)
                
        return coordinates
    
    def plotEventDisplay(self,data, color_map=plt.cm.plasma, color_norm=colors.LogNorm(), color_label=None,fig_width=None, show_zero=False, style=None ):
        #simplified version of WatChMaL plotting borrowing at lot of options and formatting
        coordinates = self.coordinates_eachChannel()
        
        #set the figure size
        if fig_width is None:
            fig_width = matplotlib.rcParams['figure.figsize'][0]
        scale = fig_width/20
        axis_ranges = np.ptp(coordinates, axis=0)
        figsize = (20*scale, 16*scale*axis_ranges[1]/axis_ranges[0])

        if not show_zero:
            data[data == 0] = np.nan
            
        if self.mask_list is not None:
            for mPMT_mask in self.mask_list:
                data[np.arange(self.nChannels)//19 == mPMT_mask] = np.nan    
            
        #set the style for plotting
        color_map = copy.copy(color_map)
        if style == "dark_background":
            edge_color = '0.35'
            color_map.set_bad(color='black')
        else:
            edge_color = '0.85'
            color_map.set_bad(color='white')
            
        if style:
            plt.style.use(style)
        
        #set the figure size and size of the PMTs in the display
        fig, ax = plt.subplots(figsize=figsize)
        fig_width = matplotlib.rcParams['figure.figsize'][0]
        # scale = fig_width/20
        
        
        #draw circles around each mPMT
        pmt_circles = [Circle((pos[0], pos[1]), radius=0.48) for pos in self.mPMT_2D_projection[:,1:3]]
        ax.add_collection(PatchCollection(pmt_circles, facecolor='none', linewidths=1*scale, edgecolors=edge_color))

        print(coordinates[:, 0].shape)
        pmts = ax.scatter(coordinates[:, 0], coordinates[:, 1], c=data, s=10*scale*scale, cmap=color_map, norm=color_norm)
        fig.colorbar(pmts, ax=ax, pad=0, label=color_label)
    
    def label_mPMT_plot(self, mPMT_label):
        
        coordinates_x = self.mPMT_2D_projection[:,1] 
        coordinates_y = self.mPMT_2D_projection[:,2]
        
        if(len(mPMT_label)!= self.nmPMTs):
            raise Exception ("label_mPMT_plot designed to have a label for each mPMT")
        
        for (x, y, label) in zip(coordinates_x,coordinates_y, mPMT_label):
            if(label is not None):
                plt.text(x, y, str(label), fontsize=10, ha='center', va='center')


    def label_mPMTs(self, mPMT_label):
        #lower left corner of each mPMT
        coordinates_x = self.mPMT_2D_projection[:,1]-0.5 
        coordinates_y = self.mPMT_2D_projection[:,2]-0.5
        
        if(len(mPMT_label)!= self.nmPMTs):
            raise Exception ("label_mPMT_plot designed to have a label for each mPMT")
        
        for (x, y, label) in zip(coordinates_x,coordinates_y, mPMT_label):
            if(label is not None):
                plt.text(x, y, str(label), fontsize=6, ha='center', va='center')
        
