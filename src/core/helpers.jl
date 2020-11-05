nw_ids(wm::AbstractWaterModel) = _IM.nw_ids(wm, :wd)
nws(wm::AbstractWaterModel) = _IM.nws(wm, :wd)

ids(wm::AbstractWaterModel, nw::Int, key::Symbol) = _IM.ids(wm, :wd, nw, key)
ids(wm::AbstractWaterModel, key::Symbol; nw::Int=wm.cnw) = _IM.ids(wm, :wd, key; nw = nw)

ref(wm::AbstractWaterModel, nw::Int = wm.cnw) = _IM.ref(wm, :wd, nw)
ref(wm::AbstractWaterModel, nw::Int, key::Symbol) = _IM.ref(wm, :wd, nw, key)
ref(wm::AbstractWaterModel, nw::Int, key::Symbol, idx) = _IM.ref(wm, :wd, nw, key, idx)
ref(wm::AbstractWaterModel, nw::Int, key::Symbol, idx, param::String) = _IM.ref(wm, :wd, nw, key, idx, param)
ref(wm::AbstractWaterModel, key::Symbol; nw::Int = wm.cnw) = _IM.ref(wm, :wd, key; nw = nw)
ref(wm::AbstractWaterModel, key::Symbol, idx; nw::Int = wm.cnw) = _IM.ref(wm, :wd, key, idx; nw = nw)
ref(wm::AbstractWaterModel, key::Symbol, idx, param::String; nw::Int = wm.cnw) = _IM.ref(wm, :wd, key, idx, param; nw = nw)

var(wm::AbstractWaterModel, nw::Int = wm.cnw) = _IM.var(wm, :wd, nw)
var(wm::AbstractWaterModel, nw::Int, key::Symbol) = _IM.var(wm, :wd, nw, key)
var(wm::AbstractWaterModel, nw::Int, key::Symbol, idx) = _IM.var(wm, :wd, nw, key, idx)
var(wm::AbstractWaterModel, key::Symbol; nw::Int = wm.cnw) = _IM.var(wm, :wd, key; nw = nw)
var(wm::AbstractWaterModel, key::Symbol, idx; nw::Int = wm.cnw) = _IM.var(wm, :wd, key, idx; nw = nw)

con(wm::AbstractWaterModel, nw::Int = wm.cnw) = _IM.con(wm, :wd; nw = nw)
con(wm::AbstractWaterModel, nw::Int, key::Symbol) = _IM.con(wm, :wd, nw, key)
con(wm::AbstractWaterModel, nw::Int, key::Symbol, idx) = _IM.con(wm, :wd, nw, key, idx)
con(wm::AbstractWaterModel, key::Symbol; nw::Int = wm.cnw) = _IM.con(wm, :wd, key; nw = nw)
con(wm::AbstractWaterModel, key::Symbol, idx; nw::Int = wm.cnw) = _IM.con(wm, :wd, key, idx; nw = nw)

sol(wm::AbstractWaterModel, nw::Int = wm.cnw) = _IM.sol(wm, :wd; nw = nw)
sol(wm::AbstractWaterModel, nw::Int, key::Symbol) = _IM.sol(wm, :wd, nw, key)
sol(wm::AbstractWaterModel, nw::Int, key::Symbol, idx) = _IM.sol(wm, :wd, nw, key, idx)
sol(wm::AbstractWaterModel, key::Symbol; nw::Int = wm.cnw) = _IM.sol(wm, :wd, key; nw = nw)
sol(wm::AbstractWaterModel, key::Symbol, idx; nw::Int = wm.cnw) = _IM.sol(wm, :wd, key, idx; nw = nw)

ismultinetwork(wm::AbstractWaterModel) = _IM.ismultinetwork(wm, :wd)
