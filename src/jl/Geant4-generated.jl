module Geant4


import Base.getindex
import Base.setindex!

using CxxWrap
import Libdl
@wrapmodule(()->"$(@__DIR__)/../deps/libjlGeant4." * Libdl.dlext)

function __init__()
    @initcxx
end

end #module
