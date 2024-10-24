import CairoMakie, GLMakie, Makie
import Statistics

function meanIntensityPlotAndIsosurface(data::AbstractArray{T,3}; func::Function=mean, fig=CairoMakie.Figure(size=(1150,300))) where T
    CairoMakie.empty!(fig) # ensure figure is empty
    imSize = size(data)
    xLabelVec = ["x", "x", "y"]
    yLabelVec = ["y", "z", "z"]
    dimVec = [3,2,1]
    sizeVec = [imSize[1:2], imSize[[1,3]], imSize[2:3]]
    axConfig = Dict(:aspect=>1, :ylabelrotation=>0, :xlabelsize=>20, :ylabelsize=>20)
    for i=1:3
        ax = CairoMakie.Axis(fig[1,i]; xlabel=xLabelVec[i], ylabel=yLabelVec[i], axConfig...)
        CairoMakie.heatmap!(ax, reshape(func(data[:,:,:], dims=dimVec[i]), sizeVec[i]), colormap=:inferno)
        CairoMakie.hidedecorations!(ax, label=false)
    end
    fgl = GLMakie.Figure(size=(300,300))
    axgl = GLMakie.Axis3(fgl[1,1], aspect=:data, azimuth=-pi/2-0.2)
    GLMakie.volume!(axgl, data, algorithm=:iso, isorange=0.15, isovalue=0.3, colormap=:viridis, colorrange=[0.0,0.2])
    GLMakie.hidedecorations!(axgl, label=false)
    GLMakie.hidespines!(axgl)
    Makie.update_state_before_display!(fgl)
    cb = Makie.colorbuffer(fgl.scene, backend=GLMakie, scalefactor=5)
    axleft, _ = CairoMakie.image(fig[1,4], cb[300:end-300,300:end-300]', axis=(yreversed=true, aspect=1))
    CairoMakie.hidedecorations!(axleft)
    CairoMakie.hidespines!(axleft)
    return fig
end