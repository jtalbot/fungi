
printf = function(s,...)
    local r = io.write(s:format(...))
    io.flush()
    return r
end

local load0 = time()
    translate(-10,5,0)
    --rotatey(math.pi/10)
    --rotatez(-math.pi/10)
    push()
        obj("scenes/xyzrgb_dragon.obj")
        translate(0,-40.5,0)
        scale(110,1,110)
        translate(-1,-2,-0.7)
            obj("scenes/cube.obj")
    pop()
    push()
        light()
        translate(0,60,0)
        scale(20,1,20)
        translate(-1,-2,-1)
            obj("scenes/cube.obj")
    pop()
printf("Load time: %.1f ms\n", (time()-load0)/1000)

local rays = 100000000

local render0 = time()
    local status = function(progress)
        printf("\r%.2f%%", progress/rays*100);
        return true
    end
    printf("\n")
    render(rays, "output/out2.exr", status)
    printf("\n")
printf("Render time: %.1f ms\n", (time()-render0)/1000)

