using TIPSBP1
using Plots
# versioninfo()

##Plot in julia - manual demo
dimer = BridgedDimer(d=1322.,x=59.,β=111.) #build dimer, can similarly substitute other default variables
test = TestData() #load TIPS-BP1' data

#Figure 2a
S1 = quintetspectrum(dimer,Spectrum(basis=:Adiabatic),Experiment(distribution="./test/testdata/powder_distribution.csv",field=test.x));
plot(test.x,test.prompt);
plot!(test.x,S1*0.8,label="Prompt")

#Figure 2b/2c
S2 = quintetspectrum(dimer,Spectrum(basis=:Diabatic,nOutput=2),Experiment(distribution="./test/testdata/powder_distribution.csv",field=test.x));
plot(test.x,test.hankel,label="Hankel");
plot!(test.x,S2[:,2]*0.8,label="⁵TT₀↔⁵TT₁")

plot(test.x,test.residual,label="Residual");
plot!(test.x,S2[:,1]*0.8,label="⁵TT₁↔⁵TT₂")

#Figure 3a
fig1, fig2, fig3 = weberplot(dimer,Spectrum(basis=:Diabatic,nOutput=2,multiOut=true),Experiment(θ=[0.,90.,90.],ϕ=[0.,0.,90.],field=test.x))
display(fig1)
display(fig2)
display(fig3)

#Figure 3d
fig = polarplot(dimer,Polarization());
display(fig)

#Figure 3e
fig = polarplot(dimer,Populations(M=0));
display(fig)
fig = polarplot(dimer,Populations(M=+1));
display(fig)
fig = polarplot(dimer,Populations(M=+2));
display(fig)