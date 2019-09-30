mutable struct Bus
   nodeID::Int
   kind::Symbol
   Pd::Float64
   Qd::Float64
   gr::Float64                  # grounding resistance
   ga::Float64                  # inverse of grounding resistance
   Pg::Float64
   Qg::Float64
   Pgmax::Float64
   Qgmax::Float64
   Pgmin::Float64
   Qgmin::Float64
   pi2::Float64                 # Objective coefficient
   pi1::Float64                 # Objective coefficient
   qobjcoeff::Float64
   Pmgcost::Float64
   Vmax::Float64
   Vmin::Float64
   Ji::Float64                  # DC induced by GIC voltagex
   coord::Vector{Float64}
   genids::Vector{Int}           # Bus can contain more than one generator
   farmids::Vector{Int}          # Bus can contain more than one farm
   outlist::Vector{Int}         # outgoing line indices
   inlist::Vector{Int}          # incoming line indices
   Gs::Float64     #shunt conductance
   Bs::Float64     #shunt susceptance
   Vm::Float64     #voltage magnitude (for PV buses)
   function Bus(nodeID, kind, Pd, Qd, Vmax, Vmin, Gs, Bs, Vm)
      if kind == 1
          k = :PQ
      elseif kind == 2
          k = :PV
      elseif kind == 3
          k = :Ref
      else
          error()
      end
      b = new(nodeID, k, Pd, Qd)
      b.gr = 0
      b.ga = 0
      b.Pg = 0
      b.Qg = 0
      b.Pgmax = 0
      b.Qgmax = 0
      b.pi1 = 0
      b.pi2 = 0
      b.Pmgcost = 0.0
      b.Vmax = Vmax
      b.Vmin = Vmin
      b.Ji = 0.0
      b.coord=[0.0, 0.0]
      b.genids = Int[]
      b.farmids = Int[]
      b.outlist = Int[]
      b.inlist = Int[]
      b.Gs = Gs
      b.Bs = Bs
      b.Vm = Vm
      return b
   end
end

# Functions for bus
function setg(b::Bus, genidx, Pg, Qg, Pgmax, Pgmin, Qgmax, Qgmin)
   b.Pg += Pg
   b.Qg += Qg
   b.Pgmax += Pgmax
   b.Pgmin += Pgmin
   b.Qgmax += Qgmax
   b.Qgmin += Qgmin
   if b.kind == :PQ
      warn("Generator $genidx was assigned to bus $(b.nodeID), but this bus has type PV")
   end
   push!(b.genids,genidx)
end

mutable struct Generator
   genID::Int
   busidx::Int
   Pg::Float64
   Qg::Float64
   Pgmax::Float64
   Pgmin::Float64
   Qgmax::Float64
   Qgmin::Float64
   pi1::Float64
   pi2::Float64
   pi3::Float64
   function Generator(genID, busidx, Pg, Qg, Pgmax, Pgmin, Qgmax, Qgmin)
      g = new(genID, busidx, Pg, Qg, Pgmax, Pgmin, Qgmax, Qgmin)
      g.pi1 = 1.0
      g.pi2 = 1.0
      g.pi3 = 1.0
      return g
   end
end

mutable struct Line
   arcID::Int
   tail::Int # the "to" node
   head::Int # the "from" node
   r::Float64 # the resistance value
   x::Float64 # the reactance value
   γ::Float64 # Conductance
   β::Float64 # Susceptance
   u::Float64 # the capacity of the line
   ratio::Float64 # turns ratio of transformer on the line
   distance_scale::Float64 # this will be used to scale u
   Imax::Float64 # max Iᵃ
   Imin::Float64 # min Iᵃ
   a::Float64 # inverse of line resistance
   Jij::Float64 # induced DC flow by GMDs
   Iij::Float64 # current causing magnetic saturation of transformer
   b_charge:: Float64 # charging susceptance
   function Line(arcID, tail, head, r, x, u, turns, d, Imax, Imin,b_charge)
      line = new(arcID, tail, head, r, x)
      line.γ = r/(r^2+x^2)
      line.β = -x/(r^2+x^2)
      line.u = u
      line.ratio = turns
      line.distance_scale = d
      line.Imax = Imax
      line.Imin = Imin
      line.a= 1/r
      line.Jij= 0
      line.Iij= Imax
      line.b_charge = b_charge
      return line
   end
end

getThermalCapacity(l::Line, mvaBase) = l.u#/mvaBase  # line limits
getSyncCapacity(l::Line, mvaBase) = l.y

mutable struct TransLine
   translineID::Int
   arcID::Int
   tail::Int # the "to" node
   head::Int # the "from" node
   r::Float64 # the resistance value
   x::Float64 # the reactance value
   γ::Float64 # Conductance
   β::Float64 # Susceptance
   u::Float64 # the capacity of the line
   ratio::Float64 # turns ratio of the line
   distance_scale::Float64 # this will be used to scale u
   Imax::Float64 # max Iᵃ
   Imin::Float64 # min Iᵃ
   a::Float64 # inverse of line resistance
   Jij::Float64 # induced DC flow by GMDs
   Iij::Float64 # current causing magnetic saturation of transformer
   b_charge::Float64
   function TransLine(translineID, arcID, tail, head, r, x, u, turns, d, Imax, Imin, b_charge)
      transline = new(translineID, arcID, tail, head, r, x)
      transline.γ = r/(r^2+x^2)
      transline.β = -x/(r^2+x^2)
      transline.u = u
      transline.ratio = turns
      transline.distance_scale = d
      transline.Imax = Imax
      transline.Imin = Imin
      transline.a= 1/r
      transline.Jij= 0.0
      transline.Iij= Imax
      transline.b_charge=b_charge
      return transline
   end
end

mutable struct Scenario
   lineIDs::Vector{Int}
   function Scenario(lineIDs)
      s = new(lineIDs)
      return s
   end
end

function GetlineID(lines, head, tail)
   lineIDs = Int[]
   for i=1:length(lines)
      if((lines[i].head == head && lines[i].tail == tail) || (lines[i].head == tail && lines[i].tail == head))
         push!(lineIDs,lines[i].arcID)  #Push all the repeated lines between given nodes.
      end
   end
   return lineIDs
end

mutable struct Farm
    μ::Float64
    σ::Float64
    bus::Int
end
