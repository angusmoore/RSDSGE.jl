function testupdateparameters(model,p,v)
    updateparameters!(model,p,v)
    return true
end

@test testupdateparameters(twoEQNK,"rho",0.95)
@test testupdateparameters(twoEQNK,["rho","phi"],[0.9,(0.8,0.6)])
