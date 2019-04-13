
printf = function(s,...)
    local r = io.write(s:format(...))
    io.flush()
    return r
end

local load0 = time()
    parsePBRT('../pbrt/pbrt-v3-scenes/barcelona-pavilion/pavilion-day.pbrt')
printf("Load time: %.1f ms\n", (time()-load0)/1000)

local rays = 800*425*16

local render0 = time()
    local status = function(progress)
        printf("\r%.2f%%", progress/rays*100);
        return true
    end
    printf("\n")
    render(rays, "output/pav.exr", status)
    printf("\n")
printf("Render time: %.1f ms\n", (time()-render0)/1000)

