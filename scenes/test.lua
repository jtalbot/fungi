
printf = function(s,...)
    return io.write(s:format(...))
end

local load0 = time()
    translate(-10,5,0)
    rotatey(math.pi/10)
    rotatez(-math.pi/10)
        obj("scenes/xyzrgb_dragon.obj")
        translate(0,-40.5,0)
        scale(110,1,110)
        translate(-1,-2,-1)
            obj("scenes/cube.obj")
printf("Load time: %.1f ms\n", (time()-load0)/1000)

local render0 = time()
    render("output/out2.exr")
printf("Render time: %.1f ms\n", (time()-render0)/1000)

