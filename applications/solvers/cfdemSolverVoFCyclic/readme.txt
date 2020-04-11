Following steps should be done to make this solver pass the compilation. 
If ESI version was used, Step 1. could be ignored.

1. 
replace the files with incompressibleTwoPhaseMixture.C/H in $(FOAM_SRC)/transportModels/incompressible/incompressibleTwoPhaseMixture/

the aim is to 
add  
Foam::tmp<Foam::scalarField>
 Foam::incompressibleTwoPhaseMixture::mu(const label patchI) const
 {
 
     return mu()().boundaryField()[patchI];
 }  
to $(FOAM_SRC)/transportModels/incompressible/
incompressibleTwoPhaseMixture/incompressibleTwoPhaseMixture.C 

add 
//- Return the dynamic laminar viscosity on patch
tmp<scalarField> mu(const label patchI) const;

to $(FOAM_SRC)/transportModels/incompressible/
incompressibleTwoPhaseMixture/incompressibleTwoPhaseMixture.H 

2. recompile the incompressibleTwoPhaseMixture

3. compile the CompressibleTwoPhaseMixtureTurbulenceModels

4. compile the expVoFDEM
