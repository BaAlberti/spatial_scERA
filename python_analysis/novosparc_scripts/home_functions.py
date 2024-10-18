import matplotlib
import plotly.graph_objs as go
import numpy as np
import pandas as pd

def embedding_3D_plotly(dataset,color,threshold=0.25,title=None,title_fig=None,screenshot=False,PCR_mode=False,pal="viridis"):
    
    """
    This function is similar to the novosparc original cluster. The aim is to plot in 3D the reconstruction for a specific gene/enhancer.
    "dataset", "color", and "title" parameters are the one already provided in the original function 
    Main changes are : 
        Changing matplotlib to plotly
        "threshold" parameter to plot expression/ probability only if it is higher than the threshold at a specific position
        "title_fig" parameter is used to specify the name of the pictures to generates.
        "screenshot" parameter is a mode telling the function if it should plot the reconstruction on the screen or take and save pictures of the reconstruction
        "PCR_mode" parameter is a mode switching how values are handled by the function depending on if you are working with expression level (False) or presence probability (True).
        "pal" parameter is to change the color scale to use for visualisation.
    """
    title = color if title is None else title
    ncolor = len(color)   
    if title_fig==None:
        title_fig=title

    # as values has to be scaled between 0 and 1, we set the normalised function of matplotlib to achieve this on our data.
    norm=matplotlib.colors.Normalize() 

    # collecting the shape of the reconstruction in 3D
    xyz = dataset.obsm['spatial']
    x = xyz[:, 0]
    y = xyz[:, 1] if xyz.shape[1] > 1 else np.ones_like(x)
    z = xyz[:, 2] if xyz.shape[1] > 1 else np.ones_like(x)

    # making a vector of normalised expression for a specific set of gene 
    for i, g in enumerate(color):
        if g in dataset.var_names:
            values = dataset[:, g].X
            values=norm(values)
            if PCR_mode==False:
                values= [ item for elem in values for item in elem]
        elif g in dataset.obs.columns:
            values = dataset.obs[g]
            values=norm(values)
            if PCR_mode==False:
                values= [ item for elem in values for item in elem]
        else:
            continue
        positive_xyz=[]
        positive_values=[]
        j=0

        # If a threshold is added, consider only the position where the expression is above this threshold
        for value in values:
            if value >= threshold:
                positive_xyz.append(xyz[j,])
                positive_values.append(value)
            j=j+1

        # Generate a plot with the shape of the embryo + the positions where the expression is above the threshold.
        fig=go.Figure(data=[go.Scatter3d(x=x,y=y,z=z,mode="markers",marker=dict(size=3.5,color='rgb(136,204,238)',opacity=0.5),name="neg cells")])
        fig.add_trace(go.Scatter3d(x=[item[0] for item in positive_xyz],y=[item[1] for item in positive_xyz],z=[item[2] for item in positive_xyz]
                                   ,mode="markers", marker_showscale=False ,marker=dict(size=4,color=np.array(positive_values),colorscale=pal,reversescale=False,cmin=0,cmax=1)
                                   ,name="pos cells"))

        # If screenshot parameter is set to True, takes pictures of three angles of the reconsruction and save them in the images folder.
        if screenshot:
            fig.update_layout(scene=dict(xaxis_title='anterior-posterior',zaxis_title='dorsal-ventral', xaxis=dict(autorange="reversed",showgrid=False,backgroundcolor="rgba(0,0,0,0)",gridcolor="white",showbackground=True,zerolinecolor="white",showticklabels=False,visible=False), yaxis=dict(showgrid=False,backgroundcolor="rgba(0,0,0,0)",gridcolor="white",showbackground=True,zerolinecolor="white",showticklabels=False,visible=False), zaxis=dict(showgrid=False,backgroundcolor="rgba(0,0,0,0)",gridcolor="white",showbackground=True,zerolinecolor="white",showticklabels=False,visible=False)), 
                              title="",
                              legend=dict(yanchor="top",y=0.99,xanchor="left",x=0.2),
                              coloraxis_colorbar=dict(title="Gene expression level", yanchor="top", y=0.8,xanchor="right",x=0.7),
                              scene_camera=dict(eye=dict(x=0, y=2, z=0)),showlegend=False)
            fig.write_image("images/"+title_fig+"_lateral_view.png",width=1000,height=500)
        
            fig.update_layout(scene=dict(xaxis_title='anterior-posterior',zaxis_title='dorsal-ventral', xaxis=dict(autorange="reversed",showgrid=False,backgroundcolor="rgba(0,0,0,0)",gridcolor="white",showbackground=True,zerolinecolor="white",showticklabels=False,visible=False), yaxis=dict(showgrid=False,backgroundcolor="rgba(0,0,0,0)",gridcolor="white",showbackground=True,zerolinecolor="white",showticklabels=False,visible=False), zaxis=dict(showgrid=False,backgroundcolor="rgba(0,0,0,0)",gridcolor="white",showbackground=True,zerolinecolor="white",showticklabels=False,visible=False)),
                              title="",
                              legend=dict(yanchor="top",y=0.99,xanchor="left",x=0.2),
                              coloraxis_colorbar=dict(title="Gene expression level", yanchor="top", y=0.8,xanchor="right",x=0.7),
                              scene_camera=dict(eye=dict(x=0, y=1, z=1.5)),showlegend=False)
            fig.write_image("images/"+title_fig+"_dorsal_view.png",width=1000,height=500)
            
            fig.update_layout(scene=dict(xaxis_title='anterior-posterior',zaxis_title='dorsal-ventral', xaxis=dict(autorange="reversed",showgrid=False,backgroundcolor="rgba(0,0,0,0)",gridcolor="white",showbackground=True,zerolinecolor="white",showticklabels=False,visible=False), yaxis=dict(showgrid=False,backgroundcolor="rgba(0,0,0,0)",gridcolor="white",showbackground=True,zerolinecolor="white",showticklabels=False,visible=False), zaxis=dict(showgrid=False,backgroundcolor="rgba(0,0,0,0)",gridcolor="white",showbackground=True,zerolinecolor="white",showticklabels=False,visible=False)),
                              title="",
                              legend=dict(yanchor="top",y=0.99,xanchor="left",x=0.2),
                              coloraxis_colorbar=dict(title="Gene expression level", yanchor="top", y=0.8,xanchor="right",x=0.7),
                              scene_camera=dict(eye=dict(x=0, y=1.5, z=-1.5)),showlegend=False)
            fig.write_image("images/"+title_fig+"_ventral_view.png",width=1000,height=500)

        else:
            fig.update_layout(scene=dict(xaxis_title='anterior-posterior',zaxis_title='dorsal-ventral', xaxis=dict(autorange="reversed",showgrid=False,backgroundcolor="rgba(0,0,0,0)",gridcolor="white",showbackground=True,zerolinecolor="white",showticklabels=False,visible=False), yaxis=dict(showgrid=False,backgroundcolor="rgba(0,0,0,0)",gridcolor="white",showbackground=True,zerolinecolor="white",showticklabels=False,visible=False), zaxis=dict(showgrid=False,backgroundcolor="rgba(0,0,0,0)",gridcolor="white",showbackground=True,zerolinecolor="white",showticklabels=False,visible=False)),
                              title=title[i],title_font=dict(family="Droid Sans",size=35),title_x=0.5,
                              scene_camera=dict(eye=dict(x=0, y=2.5, z=0)),
                              height=400, showlegend=False)
            fig.show()    




def find_index(dataset,list_cell):
    """
    This function extract the position index of specified cells in a dataset.
    """
    query_index=[]
    dataset_cells=list(dataset.obs_names)
    for cell_name in list_cell:
        ind=dataset_cells.index(cell_name)
        query_index.append(ind)
    return query_index

def prob_distributions(tissue,reconstruction,index_list,elem="sum"):
    """
    This function sum the presence probability of all the positions extracted by the find_index function and add this column in the reconstruction.
    "issue" parameter is the tissue generated after the reconstruction (not the dataset but the "tissue" object)
    "reconstruction" parameter is the dataset after reconstruction
    "index_list" parameter is the list of index generated by the find_index function
    "elem" parameter is the name of the column added to the dataset
    """
    gw = tissue.gw
    ngw = (gw.T / gw.sum(1)).T
    cell_idx = index_list
    cell_prb_cols = ['cell %d' % i for i in cell_idx]
    individual_columns = pd.DataFrame(ngw.T[:, cell_idx], columns=cell_prb_cols)
    sum_column=pd.DataFrame(individual_columns.sum(axis=1),columns=[elem]) #Somme par position
    reconstruction.obs=sum_column
    return reconstruction