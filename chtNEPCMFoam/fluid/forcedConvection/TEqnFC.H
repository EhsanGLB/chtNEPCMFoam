
{
    volScalarField kappaEff("kappaEff", (turb.nut()/Prt) * rho * Cp + kappa);
    kappaEff.correctBoundaryConditions();

    Kappa = kappaEff;

    fvScalarMatrix TEqn
    (
	fvm::ddt(rho*Cp, T)
      + rho*Cp*fvm::div(phi, T)
      - fvm::laplacian(kappaEff, T)
    );

    TEqn.relax();
    TEqn.solve();

}
