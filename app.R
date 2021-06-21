library(dplyr)
library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(rgl)
library(readobj)
library(RColorBrewer)
library(plotly)
library(geomorph)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
rm(list = ls())
setwd("~/project/st/data/E16-18_a_fig_0619/")
#read meta file

df = read.csv("merge_slice.obs.csv",sep = ",",stringsAsFactors = F)
# read marker from h5ad 
#Convert("merge_slice.h5ad",dest = "h5seurat", overwrite = TRUE)
seu = LoadH5Seurat("merge_slice.h5seurat")  
data = seu@assays$RNA@data
meta = seu@meta.data
rm(seu)

data =  as.data.frame(t(data))
meta$coord = rownames(meta)
data$coord = rownames(data)
#merge
df = dplyr::left_join(meta,data,by = "coord")
rm(meta)
rm(data)
##normalize coord in df
df$x = (df$new_x-min(df$new_x))/20
df$y = (df$new_y-min(df$new_y))/20
df$z = df$new_z/10000


#model
mesh_lst = read.obj("organ.obj")[[1]]
#get shell
shell = read.obj("shell.obj")[[1]][[1]]
###normalize the coord
xmin = min(shell$positions[3,])
ymin = min(shell$positions[2,])
zmin = min(shell$positions[1,])
shell =  list(vt = data.frame(x = shell$positions[3,]-xmin,
                              y = shell$positions[2,]-ymin,
                              z = shell$positions[1,]-zmin),
              i = shell$indices[1,], j = shell$indices[2,], k = shell$indices[3,])
##get organ
df_mesh_1st = lapply(mesh_lst, function(m){
  list(vt = data.frame(x = m$positions[3,]-xmin,
                       y = m$positions[2,]-ymin,
                       z = m$positions[1,]-zmin),
       i = m$indices[1,], j = m$indices[2,], k = m$indices[3,]
  )
})
temp = length(df_mesh_1st)
cols =  paste0("cluster",c(1:temp))
cols = c("shell",cols)
## set scene
scene = list(yaxis = list(title = '',range = round(range(shell$vt$y),0)),
             zaxis = list(title = '',range = round(range(shell$vt$z),0), scaleanchor = "y"),
             xaxis = list(title = '',range = round(range(shell$vt$x),0), scaleanchor = "y"),
             aspectratio = list(x=(round(range(shell$vt$x),0)[2]/ round(range(shell$vt$y),0)[2]), 
                                y=1, 
                                z=(round(range(shell$vt$z),0)[2]/ round(range(shell$vt$y),0)[2])))
col.base <- colorRampPalette(brewer.pal(8, "Set1"))(18)
app <- Dash$new()
i =1
df$mark = as.numeric(df$`128up`)
available_cluster = unique(as.numeric(df$clusters_0.3))
option_cluster <- lapply(available_cluster,
                         function(available_cluster) {
                           list(label = available_cluster,
                                value = available_cluster)
                         }
)


defaultPlot <- plot_ly()%>%
  add_markers( type = "scatter3d", 
               mode = "markers", 
               data = df[df$mark>3,c("x","y","z","mark")], 
               x = ~x, 
               y = ~y, 
               z = ~z,
               marker = list(opacity=0.6, 
                             size=1, color=~mark, 
                             colorscale ="Bluered")
  )%>%
  add_trace(type = 'mesh3d',
            data = shell$vt,
            x = ~x,
            y = ~y,
            z = ~z,
            i = shell$i, j = shell$j, k = shell$k,
            facecolor = rep("gray",length(shell$i)),
            opacity = 0.1,
            name = cols[1]
  )%>%  add_trace(type = 'mesh3d',
                  data = df_mesh_1st[[i]]$vt,
                  x = ~x,
                  y = ~y,
                  z = ~z,
                  i = df_mesh_1st[[i]]$i, j = df_mesh_1st[[i]]$j, k = df_mesh_1st[[i]]$k,
                  facecolor = rep(col.base[i],length(df_mesh_1st[[i]]$i)),
                  opacity = 1,name = cols[i+2]
  )
app$layout(
  htmlDiv(
    list(
      htmlH5("Input your gene:"),
      dccInput(id = "graphTitle",
               value = "128up",
               type = "text"),
      htmlDiv(id = "outputID"),
      htmlH5("Do you want to show cluster ?"),
      dccDropdown(
        id = 'cluster',
        options = option_cluster,
        value = 1
      ),
      dccGraph(id = "giraffe",
               figure = defaultPlot
      )
    )
  )
)

app$callback(output = list(id = "giraffe", property = "figure"),
             params = list(input("graphTitle", "value"),
                           input("cluster", "value")),
             function(mark,i){
               i = as.numeric(i)
               df$mark = df[,mark]
               p <- plot_ly()%>%
                 add_markers( type = "scatter3d", 
                              mode = "markers", 
                              data = df[df$mark>3,c("x","y","z","mark")], 
                              x = ~x, 
                              y = ~y, 
                              z = ~z,
                              marker = list(opacity=0.6, 
                                            size=1, color=~mark, 
                                            colorscale ="Bluered")
                 ) %>% add_trace(type = 'mesh3d',
                                 data = shell$vt,
                                 x = ~x,
                                 y = ~y,
                                 z = ~z,
                                 i = shell$i, j = shell$j, k = shell$k,
                                 facecolor = rep("gray",length(shell$i)),
                                 opacity = 0.2,
                                 name = cols[1]
                 )%>%layout(scene = scene)%>%
                 add_trace(type = 'mesh3d',
                           data = df_mesh_1st[[i]]$vt,
                           x = ~x,
                           y = ~y,
                           z = ~z,
                           i = df_mesh_1st[[i]]$i, j = df_mesh_1st[[i]]$j, k = df_mesh_1st[[i]]$k,
                           facecolor = rep(col.base[i],length(df_mesh_1st[[i]]$i)),
                           opacity = 1,name = cols[i+2]
                 )
               return(p)
             })

#df$CG5028
app$run_server(host = "0.0.0.0", port = Sys.getenv('PORT', 8050),showcase = TRUE)
