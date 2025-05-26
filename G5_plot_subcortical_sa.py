
import numpy as np
import nibabel as nib
from nilearn import datasets, plotting, image
import pandas as pd
import os
import matplotlib.pyplot as plt  # Import Matplotlib for saving

from nilearn.image import new_img_like

df_k=pd.read_csv('/15q11.2_extreme_deviations_subcortical.csv')  # Using tab as the separator)

# Load the atlas (e.g., Harvard-Oxford subcortical atlas)
atlas_path='/HarvardOxford-sub-maxprob-thr25-2mm.nii.gz'
# Load the atlas using nibabel
atlas_img = nib.load(atlas_path)
# Get the atlas data
atlas_data = atlas_img.get_fdata()
unique_labels = np.unique(atlas_data)
atlas_labels = unique_labels
# Create an empty array for storing the data to be plotted
data_to_plot = np.zeros(atlas_data.shape)

df_del=df_k[df_k['cnv'] == '15q11.2del']
df_dup=df_k[df_k['cnv'] == '15q11.2dup']
# Replace this with your actual data

roi_to_label = {
   #  "Left-Cerebellum-White-Matter":1,
   #  "Left-Cerebellum-Cortex":2,
     "left lateral ventrical":3,
     "Left-Thalamus-Proper":4, 
   	"Left-Caudate":5,
    "Left-Putamen":6,
     "Left-Pallidum":7,
    "Brainstem":8,
     'Left-Hippocampus':9, 
     "Left-Amygdala":10, 
     "Left-Accumbens-area":11, 
   #  "Right-Cerebellum-White-Matter":12,
   #  "Right-Cerebellum-Cortex":13,
     "left lateral ventrical":14,
     "Right-Thalamus-Proper":15, 
    	"Right-Caudate":16,
     "Right-Putamen":17,
     "Right-Pallidum":18,
     'Right-Hippocampus':19, 
     "Right-Amygdala":20, 
     "Right-Accumbens-area":21, 

    # Extend this dictionary based on actual documentation
}



df_plot=df_dup ###
df_plot_cnv='dup'###change
df_plot['label_number'] = df_plot['roi'].map(roi_to_label) ###change
data_values=df_plot###change

# Traverse the labels
for label_index, label_name in enumerate(atlas_labels):
    if label_name in data_values['label_number'].values:
        value = data_values.loc[data_values['label_number'] == label_name, 'p_n_outliers'].values[0]
        value=value*100
        data_to_plot[atlas_data == label_index] = value
        
        
# Create a new image with mapped data
mapped_img = new_img_like(atlas_img, data_to_plot)

# Plot the volume data as a statistical map
plotting.plot_stat_map(mapped_img, title='Extreme Neg Deviations of Subcortical of '+df_plot_cnv,
                       cmap='rocket', cut_coords=(0, 0, 0), display_mode='ortho',
                       vmin=0, 
                       vmax=8)
fig = plt.gcf()
plotting.show()
fig.savefig('+'.png', dpi=300, bbox_inches='tight')
plt.close()
                       #output_file='/subcortical_plot_deletion.png')  # Only show top 5% for clarity







######for SA##########################################
area_z=pd.read_csv('SA_Z_transfer.csv')
col_th=area_all.columns
col_th = pd.Index(col_th.tolist()[1:])


# Assign the columns back 
if len(col_th) == area_z.shape[1]:
    area_z.columns = col_th
#area_z['SA'] = area_z.sum(axis=1) 
area_te=pd.read_csv(wdir+'/te_sa_noid.csv')
area_z=area_z.reset_index()
#area_te=area_te.reset_index()

###use te_cnv
area_te=te_cnv[[ 'ID', 'Pathogenic_CNVs']]
area_z = area_z.rename(columns={col: col.replace('_area', '') for col in area_z.columns if col.endswith('_area')})


#area_z=merged_sa
df_rh = area_z[[col for col in area_z.columns if col.startswith('rh_')]]
# Create a DataFrame with columns starting with 'lh_'
df_lh = area_z[[col for col in area_z.columns if col.startswith('lh_')]]

area_rh = df_rh.rename(columns={col: col.replace('rh_', '', 1) for col in df_rh.columns if col.startswith('rh_')})
area_lh = df_lh.rename(columns={col: col.replace('lh_', '', 1) for col in df_lh.columns if col.startswith('lh_')})

area_rh=pd.merge(area_rh, area_te, left_index=True, right_index=True)
area_lh=pd.merge(area_lh, area_te, left_index=True, right_index=True)



####just cortical deviation-no direction

def run_deviations(df_no_duplicates, direction):
    roi_ids=['G_and_S_frontomargin', 'G_and_S_occipital_inf', 'G_and_S_paracentral',
           'G_and_S_subcentral', 'G_and_S_transv_frontopol', 'G_and_S_cingul-Ant',
           'G_and_S_cingul-Mid-Ant', 'G_and_S_cingul-Mid-Post',
           'G_cingul-Post-dorsal', 'G_cingul-Post-ventral', 'G_cuneus',
           'G_front_inf-Opercular', 'G_front_inf-Orbital', 'G_front_inf-Triangul',
           'G_front_middle', 'G_front_sup', 'G_Ins_lg_and_S_cent_ins',
           'G_insular_short', 'G_occipital_middle', 'G_occipital_sup',
           'G_oc-temp_lat-fusifor', 'G_oc-temp_med-Lingual',
           'G_oc-temp_med-Parahip', 'G_orbital', 'G_pariet_inf-Angular',
           'G_pariet_inf-Supramar', 'G_parietal_sup', 'G_postcentral',
           'G_precentral', 'G_precuneus', 'G_rectus', 'G_subcallosal',
           'G_temp_sup-G_T_transv', 'G_temp_sup-Lateral', 'G_temp_sup-Plan_polar',
           'G_temp_sup-Plan_tempo', 'G_temporal_inf', 'G_temporal_middle',
           'Lat_Fis-ant-Horizont', 'Lat_Fis-ant-Vertical', 'Lat_Fis-post',
           'Pole_occipital', 'Pole_temporal', 'S_calcarine', 'S_central',
           'S_cingul-Marginalis', 'S_circular_insula_ant', 'S_circular_insula_inf',
           'S_circular_insula_sup', 'S_collat_transv_ant', 'S_collat_transv_post',
           'S_front_inf', 'S_front_middle', 'S_front_sup', 'S_interm_prim-Jensen',
           'S_intrapariet_and_P_trans', 'S_oc_middle_and_Lunatus',
           'S_oc_sup_and_transversal', 'S_occipital_ant', 'S_oc-temp_lat',
           'S_oc-temp_med_and_Lingual', 'S_orbital_lateral',
           'S_orbital_med-olfact', 'S_orbital-H_Shaped', 'S_parieto_occipital',
           'S_pericallosal', 'S_postcentral', 'S_precentral-inf-part',
           'S_precentral-sup-part', 'S_suborbital', 'S_subparietal',
           'S_temporal_inf', 'S_temporal_sup', 'S_temporal_transverse',
           ]

    df_del=df_no_duplicates[df_no_duplicates['Pathogenic_CNVs'] == '']
    df_dup=df_no_duplicates[df_no_duplicates['Pathogenic_CNVs'] == '']
    
    #df_hc=df_no_duplicates[df_no_duplicates['Pathogenic_CNVs'] == 'hc']
    
    groups=[df_del, df_dup ]
    groups_s=['', '']
    df_a=[]
    df_k=[]
    #seriesmap=X1[0]
    for k, m in zip( groups, groups_s):
        for i in roi_ids:
            i_df= np.array(k[i]) #each brain region
            p_outliers=0
            n_outliers=0
            n_subject=len(k.index)
            for j in i_df: #go through each Z scores if above or below +/-1.96 
                if j >= 1.96:
                    p_outliers+=1
                elif j <= -1.96:
                    n_outliers+=1
            
            df_a.append((m, i, p_outliers, p_outliers/n_subject, n_outliers, n_outliers/n_subject))
        
    
        #print(f'number of outliers for {i}:', n_outliers, 'percentage:',n_outliers/142978 )
    df_k=pd.DataFrame(df_a, columns=('cnv','roi', 'num_p_outliers', 'p_p_outliers','num_n_outliers', 'p_n_outliers'))
    df_k.to_csv(''.csv')  # Using tab as the separator)
    return df_k
    

deviations_rh=run_deviations(area_rh, 'rh')
deviations_lh=run_deviations(area_lh, 'lh')
#deviations_rh = deviations_rh.drop(deviations_rh.index[-1])
#deviations_lh = deviations_lh.drop(deviations_lh.index[-1])



def plot_dev_sa(deviations_rh, deviations_lh, cnv):
    deviations_rh=deviations_rh[deviations_rh['cnv'] == '2q13'+cnv]
    deviations_lh=deviations_lh[deviations_lh['cnv'] == '2q13'+cnv]
    
    
    fs_home = '/'
    lh_annot = fs_home + 'lh.aparc.a2009s.annot'
    rh_annot = fs_home + 'rh.aparc.a2009s.annot'
    
    # Load the annotation files
    lh_labels, lh_ctab, lh_names = nib.freesurfer.read_annot(lh_annot)
    rh_labels, rh_ctab, rh_names = nib.freesurfer.read_annot(rh_annot)
    
    #rh_label_mapping = {labels: name.decode('utf-8') for labels, name in enumerate(rh_names)}
    
    # Load FreeSurfer fsaverage surface meshes
    fsaverage = datasets.fetch_surf_fsaverage()
    
    # Path to your right hemisphere Destrieux annotation file
    #rh_annot = '/path/to/your/rh.aparc.a2009s.annot'  # Use your local path
    
    # Read the annotation file
    #rh_labels, rh_ctab, rh_names = nib.freesurfer.read_annot(rh_annot)
    
    decoded_rh_names = [name.decode('utf-8') for name in rh_names]

    # Create a dictionary mapping from region names to labels
    name_to_label_map = {name: label for label, name in zip(np.unique(rh_labels), decoded_rh_names)}

    
    # Create a surface array corresponding to your Destrieux labels
    surface_values = np.zeros_like(rh_labels, dtype=float)
    
    # Populate the surface data array with values from your DataFrame
    for idx, row in deviations_rh.iterrows():
        region_name = row['roi']
        value = row['p_n_outliers']*100
        if region_name in name_to_label_map:
            label = name_to_label_map[region_name]
            surface_values[rh_labels == label] = value

    
    # Plot the surface with mapped values
    plotting.plot_surf_stat_map(fsaverage['pial_right'], stat_map=surface_values,
                                hemi='right',cmap='rocket', 
                                title='Right Hemisphere Extreme Deviation of '+cnv,
                                bg_map=fsaverage['sulc_right'],vmax=8 ,colorbar=True,
                                output_file='+cnv+'_.png')  # Only show top 5% for clarity

    
        # Display the plot
    plotting.show()
    #plt.savefig(os.path.join(figure_dir,cnv+'_15q11_2_right.png'), dpi=300,)
    plt.close()
    # Final adjustments and show plot
    
    
    ###left hemisphere
    #lh_label_mapping = {labels: name.decode('utf-8') for labels, name in enumerate(lh_names)}
    
    # Load FreeSurfer fsaverage surface meshes
    fsaverage = datasets.fetch_surf_fsaverage()
    
    # Path to your right hemisphere Destrieux annotation file
    #rh_annot = '/path/to/your/rh.aparc.a2009s.annot'  # Use your local path
    
    # Read the annotation file
    lh_labels, lh_ctab, lh_names = nib.freesurfer.read_annot(lh_annot)
    
    decoded_lh_names = [name.decode('utf-8') for name in lh_names]

    # Create a dictionary mapping from region names to labels
    name_to_label_map = {name: label for label, name in zip(np.unique(lh_labels), decoded_lh_names)}

    
    # Create a surface array corresponding to your Destrieux labels
    surface_values = np.zeros_like(lh_labels, dtype=float)
    
    # Populate the surface data array with values from your DataFrame
    for idx, row in deviations_lh.iterrows():
        region_name = row['roi']
        value = row['p_n_outliers']*100
        if region_name in name_to_label_map:
            label = name_to_label_map[region_name]
            surface_values[rh_labels == label] = value

    
    # Plot the surface with mapped values
    plotting.plot_surf_stat_map(fsaverage['pial_left'], stat_map=surface_values,
                                hemi='left',cmap='rocket', 
                                title='Left Hemisphere Extreme Deviation of '+cnv,
                                bg_map=fsaverage['sulc_left'],vmax=8 ,colorbar=True,
                                output_file='+cnv+'.png')  # Only show top 5% for clarity

    
        # Display the plot
    
    #plt.savefig(os.path.join(figure_dir,cnv+'_left.png'), dpi=300,)
    plotting.show()
    plt.close()
    return cnv

cnv='del' #or dup   
cnv=plot_dev_sa(deviations_rh, deviations_lh, cnv)

