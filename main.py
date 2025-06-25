import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from EventDisplay import EventDisplay
import matplotlib.colors as colors

#1) make an instance of the event display class
eventDisplay = EventDisplay() 

#2) start by loading in the CSV file for how the mPMTs are mapped to 2d event display
#unwraps based on the mPMT slot ID 
eventDisplay.load_mPMT_positions('mPMT_2D_projection_angles.csv')

#mask out mPMT slots - newer WCSim doesn't have these mPMTs loaded  
#WCTE slot numbering
# eventDisplay.mask_mPMTs([45,77,79,27,32,85,91,99,12,14,16,18])
#WCSim container numbering
# eventDisplay.mask_mPMTs([20,73,38,49,55,65,67,33,71,92,101,95])

#3) load the WCSim mapping tube no to slot number
#for WCSim using the numpy output we need the mapping between the tube_number in WCSim and the slot and mPMT number in the detector
#this can be obtained from the geofile that WCSim produces 
#This changes if the CDS is implemented or not
eventDisplay.load_wcsim_tubeno_mapping("geofile_WCTE.txt")

#4) debug by plotting some geometry files 
data = np.load("../WCSim_work_dir/wcsim_wCDS_mu-_Beam_400MeV_20cm_0000.geo.npz", allow_pickle=True)
tube_no = data["tube_no"]
position = data["position"]

axes = ["x","y","z"]

for i, axis in enumerate(axes):
    
    #map from the tube number to the mPMT slot and position number
    data_to_plot = position[:,i]
    mPMT_id, PMT_pos = eventDisplay.map_wcsim_tubeno_to_slot_pmt_id(tube_no)
    event_disp_pmt_no = (mPMT_id*19)+PMT_pos
    data_to_plot = eventDisplay.process_data(event_disp_pmt_no,data_to_plot)
    eventDisplay.plotEventDisplay(data_to_plot,color_norm=colors.Normalize(), style= "dark_background")

    eventDisplay.label_mPMTs(np.arange(0,106))
    plt.savefig("Event_display_PMT_"+axis+"_pos.png")
    plt.close()


#load the data to plot
data = np.load("../WCSim_work_dir/wcsim_wCDS_mu-_Beam_400MeV_20cm_0000.npz", allow_pickle=True)
eventID =0
tube_no = data["digi_hit_pmt"][eventID]
data_to_plot = data["digi_hit_charge"][eventID]

# tube_no = data["tube_no"]
# position = data["position"]

#map from the tube number to the mPMT slot and position number
mPMT_id, PMT_pos = eventDisplay.map_wcsim_tubeno_to_slot_pmt_id(tube_no)
event_disp_pmt_no = (mPMT_id*19)+PMT_pos
data_to_plot = eventDisplay.process_data(event_disp_pmt_no,data_to_plot)
eventDisplay.plotEventDisplay(data_to_plot,color_norm=colors.Normalize(), style= "dark_background")

# eventDisplay.label_mPMTs(np.arange(0,106))
plt.savefig("Event_display_ev"+str(eventID)+".png")
plt.close()
    

