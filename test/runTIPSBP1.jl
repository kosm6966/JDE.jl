using TIPSBP1
using Plots

##Export to text files in ./Data/ folder
TIPSBP1.figure2()
TIPSBP1.figure3a()
TIPSBP1.figure3de()


##This was the old code that generates figures
test,experiment,dimer = initializetest()

#Figure 2a
S = quintetspectrum(dimer,Spectrum(basis=:Adiabatic),experiment);
plot(test.x,test.prompt);
plot!(test.x,S*0.8,label="Prompt")

#Figure 2b/2c
S = quintetspectrum(dimer,Spectrum(nOutput=2),experiment);
plot(test.x,test.hankel,label="Hankel");
plot!(test.x,S[:,2]*0.8,label="⁵TT₀↔⁵TT₁")
plot(test.x,test.residual,label="Residual");
plot!(test.x,S[:,1]*0.8,label="⁵TT₁↔⁵TT₂")

#Figure 3a
fig1, fig2, fig3 = weberplot(dimer,Spectrum(nOutput=2,multiOut=true),Experiment(θ=[0.,90.,90.],ϕ=[0.,0.,90.],field=test.x))
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