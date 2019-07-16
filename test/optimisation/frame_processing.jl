using CellularAutomataBase, Dispersal, Test, Colors
using CellularAutomataBase: frametoimage

@testset "image processor colors regions by fit" begin

    init =  [1.0  1.0  0.5
             0.0  0.0  0.0
             0.0  0.0  0.0]

    regionlookup = [1  2  3
                    1  2  3
                    2  2  0]

    occurance = Bool[0  1  0  # 1
                     1  0  1  # 2
                     1  1  0] # 3

    image = [RGB24(1.0,0.0,0.0)  RGB24(1.0,1.0,1.0)  RGB24(0.502,0.502,0.502)
             RGB24(1.0,0.0,1.0)  RGB24(1.0,1.0,0.0)  RGB24(1.0,1.0,0.0)
             RGB24(1.0,1.0,0.0)  RGB24(1.0,1.0,0.0)  RGB24(0.0,1.0,1.0)]

    framesperstep = 2
    truepostivecolor =   (1.0, 1.0, 1.0)
    falsepositivecolor = (1.0, 0.0, 0.0)
    truenegativecolor =  (1.0, 0.0, 1.0)
    falsenegativecolor = (1.0, 1.0, 0.0)
    maskcolor =          (0.0, 1.0, 1.0)
    objective = RegionObjective(0.0, regionlookup, occurance, framesperstep)

    processor = ColorRegionFit(objective, truepostivecolor, falsepositivecolor,
                               truenegativecolor, falsenegativecolor, maskcolor)

    output = ArrayOutput(init, 1)

    @test image == frametoimage(processor, output, init, 1)
end

