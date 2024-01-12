module Geant4


import Base.getindex
import Base.setindex!

using CxxWrap
@wrapmodule(()->"libjlGeant4")

function __init__()
    @initcxx
end

end #module
