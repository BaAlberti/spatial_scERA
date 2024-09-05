def embedding_3D_plotly(dataset,color,threshold=0.25,title=None,title_fig=None,screenshot=False,PCR_mode=False,pal="viridis"):
    
    """
    Plots fields (color) of Scanpy AnnData object on spatial coordinates
    dataset -- Scanpy AnnData with 'spatial' matrix in obsm containing the spatial coordinates of the tissue
    color -- a list of fields - gene names or columns from obs to use for color
    """
    title = color if title is None else title
    ncolor = len(color)   
    if title_fig==None:
        title_fig=title

    norm=matplotlib.colors.Normalize()
    
    # remove if not necessary   - make a real normalisation between 0 and 1    
    #def norm(x):
        #xnorm=(x-np.min(x))/(np.max(x)-np.min(x))
        #return(xnorm)
    
    xyz = dataset.obsm['spatial']
    x = xyz[:, 0]
    y = xyz[:, 1] if xyz.shape[1] > 1 else np.ones_like(x)
    z = xyz[:, 2] if xyz.shape[1] > 1 else np.ones_like(x)

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
        
        for value in values:
            if value >= threshold:
                positive_xyz.append(xyz[j,])
                positive_values.append(value)
            j=j+1
            
        fig=go.Figure(data=[go.Scatter3d(x=x,y=y,z=z,mode="markers",marker=dict(size=3.5,color='rgb(136,204,238)',opacity=0.5),name="neg cells")])
        fig.add_trace(go.Scatter3d(x=[item[0] for item in positive_xyz],y=[item[1] for item in positive_xyz],z=[item[2] for item in positive_xyz]
                                   ,mode="markers", marker_showscale=False ,marker=dict(size=4,color=np.array(positive_values),colorscale=pal,reversescale=False,cmin=0,cmax=1) # plotly3/haline/thermal/cividis/magenta
                                   ,name="pos cells"))

        if screenshot:
            fig.update_layout(scene=dict(xaxis_title='anterior-posterior',zaxis_title='dorsal-ventral',
                                         xaxis=dict(autorange="reversed",showgrid=False,backgroundcolor="rgba(0,0,0,0)",gridcolor="white",showbackground=True,zerolinecolor="white",showticklabels=False,visible=False),
                                         yaxis=dict(showgrid=False,backgroundcolor="rgba(0,0,0,0)",gridcolor="white",showbackground=True,zerolinecolor="white",showticklabels=False,visible=False),
                                         zaxis=dict(showgrid=False,backgroundcolor="rgba(0,0,0,0)",gridcolor="white",showbackground=True,zerolinecolor="white",showticklabels=False,visible=False)),
                              title="",
                              legend=dict(yanchor="top",y=0.99,xanchor="left",x=0.2),
                              coloraxis_colorbar=dict(title="Gene expression level", yanchor="top", y=0.8,xanchor="right",x=0.7),
                              scene_camera=dict(eye=dict(x=0, y=2, z=0)),showlegend=False)
            fig.write_image("images/"+title_fig+"_lateral_view.png",width=1000,height=500)
            fig.update_layout(scene=dict(xaxis_title='anterior-posterior',zaxis_title='dorsal-ventral',
                                         xaxis=dict(autorange="reversed",showgrid=False,backgroundcolor="rgba(0,0,0,0)",gridcolor="white",showbackground=True,zerolinecolor="white",showticklabels=False,visible=False),
                                         yaxis=dict(showgrid=False,backgroundcolor="rgba(0,0,0,0)",gridcolor="white",showbackground=True,zerolinecolor="white",showticklabels=False,visible=False),
                                         zaxis=dict(showgrid=False,backgroundcolor="rgba(0,0,0,0)",gridcolor="white",showbackground=True,zerolinecolor="white",showticklabels=False,visible=False)),
                              title="",
                              legend=dict(yanchor="top",y=0.99,xanchor="left",x=0.2),
                              coloraxis_colorbar=dict(title="Gene expression level", yanchor="top", y=0.8,xanchor="right",x=0.7),
                              scene_camera=dict(eye=dict(x=0, y=1, z=1.5)),showlegend=False)
            fig.write_image("images/"+title_fig+"_dorsal_view.png",width=1000,height=500)
            fig.update_layout(scene=dict(xaxis_title='anterior-posterior',zaxis_title='dorsal-ventral',
                                         xaxis=dict(autorange="reversed",showgrid=False,backgroundcolor="rgba(0,0,0,0)",gridcolor="white",showbackground=True,zerolinecolor="white",showticklabels=False,visible=False),
                                         yaxis=dict(showgrid=False,backgroundcolor="rgba(0,0,0,0)",gridcolor="white",showbackground=True,zerolinecolor="white",showticklabels=False,visible=False),
                                         zaxis=dict(showgrid=False,backgroundcolor="rgba(0,0,0,0)",gridcolor="white",showbackground=True,zerolinecolor="white",showticklabels=False,visible=False)),
                              title="",
                              legend=dict(yanchor="top",y=0.99,xanchor="left",x=0.2),
                              coloraxis_colorbar=dict(title="Gene expression level", yanchor="top", y=0.8,xanchor="right",x=0.7),
                              scene_camera=dict(eye=dict(x=0, y=1.5, z=-1.5)),showlegend=False)
            fig.write_image("images/"+title_fig+"_ventral_view.png",width=1000,height=500)

        else:
            fig.update_layout(scene=dict(xaxis_title='anterior-posterior',zaxis_title='dorsal-ventral',
                                         xaxis=dict(autorange="reversed",showgrid=False,backgroundcolor="rgba(0,0,0,0)",gridcolor="white",showbackground=True,zerolinecolor="white",showticklabels=False,visible=False),
                                         yaxis=dict(showgrid=False,backgroundcolor="rgba(0,0,0,0)",gridcolor="white",showbackground=True,zerolinecolor="white",showticklabels=False,visible=False),
                                         zaxis=dict(showgrid=False,backgroundcolor="rgba(0,0,0,0)",gridcolor="white",showbackground=True,zerolinecolor="white",showticklabels=False,visible=False)),
                              title=title[i],title_font=dict(family="Droid Sans",size=35),title_x=0.5,
                              scene_camera=dict(eye=dict(x=0, y=2.5, z=0)),
                              height=400, showlegend=False)
            fig.show()    




def find_index(dataset,list_cell):
    query_index=[]
    dataset_cells=list(dataset.obs_names)
    for cell_name in list_cell:
        ind=dataset_cells.index(cell_name)
        query_index.append(ind)
    return query_index

def prob_distributions(tissue,reconstruction,index_list,elem="sum"):
    gw = tissue.gw
    ngw = (gw.T / gw.sum(1)).T
    cell_idx = index_list
    cell_prb_cols = ['cell %d' % i for i in cell_idx]
    individual_columns = pd.DataFrame(ngw.T[:, cell_idx], columns=cell_prb_cols)
    sum_column=pd.DataFrame(individual_columns.sum(axis=1),columns=[elem]) #Somme par position
    reconstruction.obs=sum_column
    return reconstruction