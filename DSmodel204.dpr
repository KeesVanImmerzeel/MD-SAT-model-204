library dsmodel204;

{.$define test} // Maakt 'testDSmodel204.log'

  {-The Wageningen Lowland Runoff Simulator (WALRUS): a lumped rainfall-runoff
    model for catchments with shallow groundwater.

    MODIFIED:
    - surface water level is fixed (not modelled);
    - Curved (instead of linear) relation between the grondwater level and
      and groundwater discharge (=topsystem 4 Triwaco).

    Ref.:

    1. Brauer, C. C., Teuling, A. J., Torfs, P. J. J. F., Uijlenhoet, R.,
    2014a. The Wageningen Lowland Runoff Simulator (WALRUS): a lumped rainfall-
    runoff model for catchments with shallow groundwater. Geosci. Model Dev.
    7, 2313–2332.
    \url{http://www.geosci-model-dev.net/7/2313/2014/gmd-7-2313-2014.pdf

    1. Brauer, C. C., Torfs, P. J. J. F., Teuling, A. J., Uijlenhoet, R.,
    2014b. The Wageningen Lowland Runoff Simulator (WALRUS): application to the
    Hupsel Brook catchment and Cabauw polder. Hydrol. Earth Syst. Sci. 18,
    4007–4028.
    \url{www.hydrol-earth-syst-sci.net/18/4007/2014/hess-18-4007-2014.pdf


    A user manual to the original R version is provided on
    www.GitHub.com/ClaudiaBrauer/WALRUS,
    containing explanations of the model structure, model variables, model
    parameters and     what these represent. It also gives detailed tips on how
    to get WALRUS running for your own catchment.
 }

  { Important note about DLL memory management: ShareMem must be the
  first unit in your library's USES clause AND your project's (select
  Project-View Source) USES clause if your DLL exports any procedures or
  functions that pass strings as parameters or function results. This
  applies to all strings passed to and from your DLL--even those that
  are nested in records and classes. ShareMem is the interface unit to
  the BORLNDMM.DLL shared memory manager, which must be deployed along
  with your DLL. To avoid using BORLNDMM.DLL, pass string information
  using PChar or ShortString parameters. }

uses
  ShareMem,
  windows, SysUtils, Classes, Math, LargeArrays,
  ExtParU, USpeedProc, uDCfunc, UdsModel, UdsModelS, xyTable, DUtils, uError,
  fmx.Dialogs;

Const
  cModelID      = 204;  {-Uniek modelnummer}

  {-Beschrijving van de array met afhankelijke variabelen}
  cNrOfDepVar   = 19;    {-Lengte van de array met afhankelijke variabelen}
  {-States}
  c_dV          = 1;     {-dV: Storage deficit (m)}
  c_dG          = 2;     {-dG: Groundwater depth (m-mv)}
  c_hQ          = 3;     {-hQ: Level quickflow reservoir = P.W (m)}
  {-Dependent variables}
  c_W           = 4;     {-W: Wetness index (-)}
  c_BdV         = 5;     {-BdV: Evapotranspiration reduction factor}
  c_dVeq        = 6;     {-dVeq: Equilibrium storage deficit (m)}
  {-External fluxes: input}
  c_P           = 7;     {-P (=PQ+PV): Precipitation (m/d)}
  c_ETpot       = 8;     {-ETpot: Potential evapotranspiration (m/d)}
  c_fXG         = 9;     {-fXG: Seepage (up/down)/extraction (m/d)}
  {-External fluxes: output}
  c_Eact        = 10;    {-Eact: Actual evapotranspiration = Beta.ETpot (m/d)}
  c_Q1          = 11;    {-Q1: Discharge of groundwater to/from primary drainage system (m/d)}
  c_Q2          = 12;    {-Q2: Discharge of groundwater to/from secundary drainage system (m/d)}
  c_Q3          = 13;    {-Q3: Discharge of groundwater to/from tertiary drainage system (m/d)}
  c_fGS         = 14;    {-fGS: Discharge of groundwater to/from drainage system (=Q1+Q2+Q3) (m/d)}
  c_fQS         = 15;    {-fQS: Quickflow (m/d)}
  {-Internal fluxes}
  c_PV          = 16;    {-PV: Precipitation into vadose zone = P.(1-W) (m/d)}
  c_PQ          = 17;    {-PQ: Precipitation into quickflow reservoir = P.W (m/d)}
  {-Other}
  c_PHIT        = 18;    {-PHIT: Phreatic head = SL-dG (m+ref)}
  c_Err         = 19;     {-Err: Water balance error (m)}

  {-Aantal keren dat een discontinuiteitsfunctie wordt aangeroepen in de procedure met
    snelheidsvergelijkingen (DerivsProc)}
  nDC = 8;

  {-Beschrijving van de array's met daarin de status van de discontinuiteitsfuncties}
  cDCfunc_Hp  = 0; {-Bij overgang van drainage naar infiltratie opeens andere weerstand}
  cDCfunc_BD1 = 1; {-Bij aansnijden slootbodemniveau primary drainage system opeens drainage}
  cDCfunc_BD2 = 2; {-Bij aansnijden slootbodemniveau secundary drainage system opeens drainage}
  cDCfunc_BD3 = 3; {-Bij aansnijden slootbodemniveau tertiary drainage system opeens drainage}
  cDCfunc_hQ  = 4; {-Als hQ >0 dan Quickflow}
  cDCfunc_dG_psi_ae = 5; {-Discontinuiteit 1 in functie dVeq_dG}
  cDCfunc_dG_0 = 6;      {-Discontinuiteit 2 in functie dVeq_dG}
  cDCfunc_W_dV_1 = 7;    {-Discontinuiteit 1 in functie W_dV}
  //cDCfunc_W_dV_2 = 8;    {-Discontinuiteit 2 in functie W_dV}

  {-Variabelen die samenhangen met het aanroepen van het model vanuit de Shell}
  cnRP    = 18;  {-Aantal RP-tijdreeksen die door de Shell moeten worden aangeleverd (in
                   de externe parameter Array EP (element EP[ indx-1 ]))}
  cnSQ    = 0;   {-Idem punt-tijdreeksen}
  cnRQ    = 0;   {-Idem lijn-tijdreeksen}

  {-Beschrijving van het eerste element van de externe parameter-array (EP[cEP0])}
  cNrXIndepTblsInEP0 = 5;    {-Aantal XIndep-tables in EP[cEP0]}
  cNrXdepTblsInEP0   = 0;    {-Aantal Xdep-tables   in EP[cEP0]}
  {-Nummering van de xIndep-tabellen in EP[cEP0]. De nummers 0&1 zijn gereserveerd}
  cTb_MinMaxValKeys   = 2;
  cTb_SoilChar        = 3;   {-Soil characteristics:
                               Column 1: b = Pore size distribution parameter (-)
                               Column 2: psi_ae = Air entry pressure (m)
                               Column 3: theta_s = Soil moisture content at saturation (-)
                               c.q. model parameters:
                               Column 4: cV = Vadose zone relaxation time (h); als groter: minder heftige respons}
  cTb_Beta            = 4; {-Table with parameters for Evapotranspiration reduction function}

  {-Beschrijving van het tweede element van de externe parameter-array (EP[cEP1])}
  {-Opmerking: table 0 van de xIndep-tabellen is gereserveerd}
  {-Nummering van de xdep-tabellen in EP[cEP1]}
  cTb_P                = 0;   {-P: Precipitation (m/d)}
  cTb_ETpot            = 1;   {-ETpot: Potential evapotranspiration (m/d)}
  cTb_fXG              = 2;   {-fXG: Seepage (up/down)/extraction (m/d)}
  cTb_Soil             = 3;   {-Soil: Soil type (-)}
  cTb_cQ               = 4;   {-cQ: Quickflow reservoir constant (d): als groter: minder heftige respons}
  cTb_SL               = 5;   {-SL: Surface level (m+ref)}
  cTb_Hp               = 6;   {-Hp: Controlled or polder water level (m+ref)}
  cTb_Wd1              = 7;   {-Wd1: Drainage resistance of primary drainage system (d)}
  cTb_Wd2              = 8;   {-Wd2: Drainage resistance of secondary drainage system (d)}
  cTb_Wd3              = 9;   {-Wd3: Drainage resistance of tertiary drainage system (d)}
  cTb_Wi1              = 10;  {-Wi1: Infiltration resistance of primary drainage system (d)}
  cTb_Wi2              = 11;  {-Wi2: Infiltration resistance of secondary drainage system (d)}
  cTb_Wi3              = 12;  {-Wi3: Infiltration resistance of tertiary drainage system (d)}
  cTb_BD1              = 13;  {-BD1: Bottom Level of primary drainage system (m+ref)         [RP10 Triwaco topsysteem 4]}
  cTb_BD2              = 14;  {-BD2: Bottom level of secondary drainage system (m+ref)       [RP11 Triwaco topsysteem 4]}
  cTb_BD3              = 15;  {-BD3: Bottom level of tertiary drainage system (m+ref)        [RP12 Triwaco topsysteem 4]}
  cTb_PHIT_init        = 16;  {-PHIT_init: Initial phreatic head (m+ref)}
  cTb_cW               = 17;  {-Wetness Index parameter (m); als groter: piet komt eerder}

  {-Model specifieke fout-codes}
//  cInvld_KeyValue1     = -999;
//  cInvld_KeyValue2     = -999;
                 {###}
  // Foutcodes model 204: -8001..-8050}
  cInvld_P         = -8001;
  cInvld_ETpot     = -8002;
  cInvld_fXG       = -8003;
  cInvld_Soil      = -8004;
  cInvld_cQ        = -8005;
  cInvld_SL        = -8006;
  cInvld_Hp        = -8007;
  cInvld_Wd1       = -8008;
  cInvld_Wd2       = -8009;
  cInvld_Wd3       = -8010;
  cInvld_Wi1       = -8011;
  cInvld_Wi2       = -8012;
  cInvld_Wi3       = -8013;
  cInvld_BD1       = -8014;
  cInvld_BD2       = -8015;
  cInvld_BD3       = -8016;
  cInvld_PHIT_init = -8017;
  cInvld_cW        = -8018;
                 {###}
//  cInvld_ParFromShell2 = -999;
                 {###}
//  cInvld_DefaultPar1   = -999;
//  cInvld_DefaultPar2   = -999;
                 {###}
//  cInvld_Init_Val1     = -999;
//  cInvld_Init_Val2     = -999;
                 {###}

var
  Indx: Integer; {-Door de Boot-procedure moet de waarde van deze index worden ingevuld,
                   zodat de snelheidsprocedure 'weet' waar (op de externe parameter-array)
				   hij zijn gegevens moet zoeken}
  ModelProfile: TModelProfile;
                 {-Object met met daarin de status van de discontinuiteitsfuncties
				   (zie nDC) }

  {-Globally defined parameters from EP[0]}

  {-Geldige range van key-/parameter/initiele waarden. De waarden van deze  variabelen moeten
    worden ingevuld door de Boot-procedure}
  cMin_Soil, cMax_Soil : Integer;
  cMin_P, cMax_P,
  cMin_ETpot, cMax_ETpot,
  cMin_fXG, cMax_fXG,
  cMin_cQ, cMax_cQ,
  cMin_SL, cMax_SL,
  cMin_Hp, cMax_Hp,
  cMin_Wd1, cMax_Wd1,
  cMin_Wd2, cMax_Wd2,
  cMin_Wd3, cMax_Wd3,
  cMin_Wi1, cMax_Wi1,
  cMin_Wi2, cMax_Wi2,
  cMin_Wi3, cMax_Wi3,
  cMin_BD1, cMax_BD1,
  cMin_BD2, cMax_BD2,
  cMin_BD3, cMax_BD3,
  cMin_cW, cMax_cW,
  cMin_PHIT_Init, cMax_PHIT_Init : Double;
  {cMin_InitVal1,  cMax_InitVal1,  ###,}

  {PHIT,}       {-phreatic head (m+ref)}
  {dG,}         {-Groundwater depth (m-mv)}
  {dV: Double;} {-Storage deficit (m)}


Procedure MyDllProc( Reason: Integer );
begin
  if Reason = DLL_PROCESS_DETACH then begin {-DLL is unloading}
    {-Cleanup code here}
	if ( nDC > 0 ) then
      ModelProfile.Free;
  end;
end;

				 
Procedure DerivsProc( var x: Double; var y, dydx: TLargeRealArray;
                      var EP: TExtParArray; var Direction: TDirection;
                      var Context: Tcontext; var aModelProfile: PModelProfile; var IErr: Integer );
{-Deze procedure verschaft de array met afgeleiden 'dydx', gegeven het tijdstip 'x' en
  de toestand die beschreven wordt door de array 'y' en de externe condities die beschreven
  worden door de 'external parameter-array EP'. Als er geen fout op is getreden bij de
  berekening van 'dydx' dan wordt in deze procedure de variabele 'IErr' gelijk gemaakt aan de
  constante 'cNoError'. Opmerking: in de array 'y' staan dus de afhankelijke variabelen,
  terwijl 'x' de onafhankelijke variabele is (meestal de tijd)}
var
  {-Sleutelwaarden (integer) afkomstig van de Shell voor de default-tabellen in EP[cEP0] }
  Soil: Integer;

  {-Waarden in default-tabellen in EP[cEP0]}
  b,            {-Pore size distribution parameter (-)}
  psi_ae,       {-Air entry pressure (mm)}
  theta_s,      {-Soil moisture content at saturation (-)}
  cW,           {-Wetness index parameter (m)}
  cV,           {-Vadose zone relaxation time (h)}
  zeta1, zeta2: {-Parameters for evapotranspiration reduction function}
  Double;
  {$ifdef test} var f: TextFile; {$endif}
  Triggered: Boolean;

  {-Parameter-waarden (Double) afkomstig van de Shell}
  P,                               {-Precipitation (m/d)}
  ETpot,                           {-Potential evapotranspiration (m/d)}
  fXG,                             {-Seepage (up/down)/extraction (m/d)}
  cQ,                              {-Quickflow reservoir constant (d)}
  SL,                              {-Surface level (m+ref)}
  Hp,                              {-controlled or polder water level (m+ref) }
  Wd1,                             {-drainage resistance of primary drainage system (d) }
  Wd2,                             {-drainage resistance of secondary drainage system (d)}
  Wd3,                             {-drainage resistance of tertiary drainage system (d)}
  Wi1,                             {-infiltration resistance of primary drainage system (d)}
  Wi2,                             {-infiltration resistance of secondary drainage system (d) }
  Wi3,                             {-infiltration resistance of tertiary drainage system (d)}
  BD1,                             {-base elevation of primary drainage system (m+ref)}
  BD2,                             {-base elevation of secondary drainage system (m+ref)}
  BD3,                             {-base elevation of tertiary drainage system (m+ref)}
  PHIT_init: Double;               {-Initial phreatic head (m+ref)}

  {-Default parameter-waarden in EP[cEP0] die alleen gevonden kunnen middels de sleutel afkomstig uit de Shell}
//  DefaultPar1, DefaultPar2,

  {-Calculated parameters}
  W,     {-Wetness index (-)}
  Beta,  {-Evapotranspiration reduction factor (-)}
  dVeq,  {-Equilibrium storage deficit (m)}
  PV,    {-Precipitation into vadose zone = P.(1-W) (m/d)}
  PQ,    {-Precipitation into quickflow reservoir = P.W (m/d)}
  Eact,  {-Actual evapotranspiration = Beta.ETpot (m)}
  fQS,   {-Quickflow (m/d)}
  {hQ,}    {-Level quickflow reservoirLevel quickflow reservoir}
  Q1,    {-Discharge of groundwater to/from primary drainage system (m/d)}
  Q2,    {-Discharge of groundwater to/from secundary drainage system (m/d)}
  Q3,    {-Discharge of groundwater to/from tertiary drainage system (m/d)}
  fGS,   {-fGS: Discharge of groundwater to/from drainage system (=Q1+Q2+Q3) (m/d)}
  dV_reken:
  Double;
  i: Integer;


//  DefaultPar2,
//  CalcPar1, CalcPar2, ### : Double;   {-Afgeleide (berekende) parameter-waarden}

Function SetParValuesFromEP0( var IErr: Integer ): Boolean;
  {-Fill globally defined parameters from EP[0]. If memory is allocated here,
    free first with 'try .. except' for those cases that the model is used repeatedly}
begin
  Result := true;
  IErr   := cNoError;
//  with EP[ cEP0 ].xInDep.Items[ cTb_Beta ] do begin
//    zeta1 := GetValue( 1, 1 );
//    zeta2 := GetValue( 1, 2 );
//  end;
//  with EP[ cEP0 ].xInDep.Items[ cTb_Beta ] do begin
//    zeta1 := GetValue( 1, 1 );
//    zeta2 := GetValue( 1, 2 );
//  end;
end;

{-Grondwaterstand (m-mv) = Surface level (SL) - freatic head (PHIT)}
Function func_dG( const SL, PHIT: Double ): Double;
begin
  Result := SL - PHIT;
end;

{-Function W(dV): Wetness index (-)}
Function func_W_dV( const dV, cW: Double ): Double;
begin
  Result := ( cos( max( min( dV, cW ), 0 )* pi / cW) / 2 + 0.5);
end;

{-Function dVeq(dG): Equilibrium storage deficit (m)}
Function func_dVeq_dG( const dG, psi_ae, theta_s, b: Double ): Double;
begin
  if( dG > psi_ae) then
    Result := ( dG
                - psi_ae /(1-b)
                - dG * Power( dG/psi_ae, -1/b )
                + psi_ae / ( 1 - b ) * Power( dG/psi_ae, (1-1/b) )
              ) * theta_s
  else if ( dG < 0) then
    Result := dG
  else
    Result := 0;
end;

{-Function fQS: Quickflow (m/d)}
Function func_fQS( const hQ, cQ: Double ): Double;
begin
  Result := hQ/cQ;
end;

{-Function beta(dV): Evapotranspiration reduction factor}
Function func_beta_dV( const dV, zeta1, zeta2: Double ): Double;
begin
  Result := ( (1-exp(-zeta1*(zeta2-dV)) ) / (1+exp(-zeta1*(zeta2-dV)) ) /2 + 0.5)
end;

Function SetKeyAndParValues( var IErr: Integer ): Boolean;
  {-Get default parameters (EP0}
  Function Get_b( const Soil: Integer ): Double;
  begin
    with EP[ cEP0 ].xInDep.Items[ cTb_SoilChar ] do
      Result := GetValue( Soil, 1 ); {row, column}
  end;

  Function Get_psi_ae( const Soil: Integer ): Double;
  begin
    with EP[ cEP0 ].xInDep.Items[ cTb_SoilChar ] do
      Result := GetValue( Soil, 2 ); {row, column}
  end;

  Function Get_theta_s( const Soil: Integer ): Double;
  begin
    with EP[ cEP0 ].xInDep.Items[ cTb_SoilChar ] do
      Result := GetValue( Soil, 3 ); {row, column}
  end;

  Function Get_cV( const Soil: Integer ): Double;
  begin
    with EP[ cEP0 ].xInDep.Items[ cTb_SoilChar ] do
      Result := GetValue( Soil, 4 ); {row, column}
  end;

  Function Get_zeta1: Double;
  begin
    with EP[ cEP0 ].xInDep.Items[ cTb_Beta ] do
      Result := GetValue( 1, 1 ); {row, column}
  end;

  Function Get_zeta2: Double;
  begin
    with EP[ cEP0 ].xInDep.Items[ cTb_Beta ] do
      Result := GetValue( 1, 2 ); {row, column}
  end;

  {-Get parameters from shell (EP1)}

  Function GetParFromShell_P( const x: Double ): Double;
  begin
    with EP[ indx-1 ].xDep do
      Result := Items[ cTb_P ].EstimateY( x, Direction );
  end;

  Function GetParFromShell_ETpot( const x: Double ): Double;
  begin
    with EP[ indx-1 ].xDep do
      Result := Items[ cTb_ETpot ].EstimateY( x, Direction );
  end;

  Function GetParFromShell_fXG( const x: Double ): Double;
  begin
    with EP[ indx-1 ].xDep do
      Result := Items[ cTb_fXG ].EstimateY( x, Direction );
  end;

  Function GetParFromShell_Soil( const x: Double ): Integer;
  begin
    with EP[ indx-1 ].xDep do
      Result := Trunc( Items[ cTb_Soil ].EstimateY( x, Direction ) );
  end;

  Function GetParFromShell_cQ( const x: Double ): Double;
  begin
    with EP[ indx-1 ].xDep do
      Result := Items[ cTb_cQ ].EstimateY( x, Direction );
  end;

  Function GetParFromShell_SL( const x: Double ): Double;
  begin
    with EP[ indx-1 ].xDep do
      Result := Items[ cTb_SL ].EstimateY( x, Direction );
  end;

  Function GetParFromShell_Hp( const x: Double ): Double;
  begin
    with EP[ indx-1 ].xDep do
      Result := Items[ cTb_Hp ].EstimateY( x, Direction );
  end;

  Function GetParFromShell_Wd1( const x: Double ): Double;
  begin
    with EP[ indx-1 ].xDep do
      Result := Items[ cTb_Wd1 ].EstimateY( x, Direction );
  end;

  Function GetParFromShell_Wd2( const x: Double ): Double;
  begin
    with EP[ indx-1 ].xDep do
      Result := Items[ cTb_Wd2 ].EstimateY( x, Direction );
  end;

  Function GetParFromShell_Wd3( const x: Double ): Double;
  begin
    with EP[ indx-1 ].xDep do
      Result := Items[ cTb_Wd3 ].EstimateY( x, Direction );
  end;

  Function GetParFromShell_Wi1( const x: Double ): Double;
  begin
    with EP[ indx-1 ].xDep do
      Result := Items[ cTb_Wi1 ].EstimateY( x, Direction );
  end;

  Function GetParFromShell_Wi2( const x: Double ): Double;
  begin
    with EP[ indx-1 ].xDep do
      Result := Items[ cTb_Wi2 ].EstimateY( x, Direction );
  end;

  Function GetParFromShell_Wi3( const x: Double ): Double;
  begin
    with EP[ indx-1 ].xDep do
      Result := Items[ cTb_Wi3 ].EstimateY( x, Direction );
  end;

  Function GetParFromShell_BD1( const x: Double ): Double;
  begin
    with EP[ indx-1 ].xDep do
      Result := Items[ cTb_BD1 ].EstimateY( x, Direction );
  end;

  Function GetParFromShell_BD2( const x: Double ): Double;
  begin
    with EP[ indx-1 ].xDep do
      Result := Items[ cTb_BD2 ].EstimateY( x, Direction );
  end;

  Function GetParFromShell_BD3( const x: Double ): Double;
  begin
    with EP[ indx-1 ].xDep do
      Result := Items[ cTb_BD3 ].EstimateY( x, Direction );
  end;

  Function GetParFromShell_cW( const x: Double ): Double;
  begin
    with EP[ indx-1 ].xDep do
      Result := Items[ cTb_cW ].EstimateY( x, Direction );
  end;


//  Function GetDefaultPar2( const KeyValue1, KeyValue2: Integer ): Double;
//  begin
//    with EP[ cEP0 ].xInDep.Items[ cTb_DefaultPar2 ] do
//      Result := GetValue( KeyValue1, KeyValue2 ); {row, column}
//  end;

  {###}

  {- User defined functions/procedures to calculate CalcPar1, CalcPar2... etc.}

  {###}


begin {-Function SetKeyAndParValues}
  Result := False;
  IErr   := cUnknownError;
  {CalcPar1 := cNoResult;}

  {###}

  Soil := GetParFromShell_Soil( x );
  if ( Soil < cMin_Soil ) or ( Soil > cMax_Soil ) then begin
    IErr := cInvld_Soil; Exit;
  end;
  b := Get_b( Soil );
  psi_ae := Get_psi_ae( Soil );
  theta_s := Get_theta_s( Soil );
  cV := Get_cV( Soil );

  zeta1 := Get_zeta1;
  zeta2 := Get_zeta2;

  {KeyValue2 := ###}
  {###}

  P := GetParFromShell_P( x );
  if ( P < cMin_P ) or ( P > cMax_P ) then begin
    IErr := cInvld_P; Exit;
  end;

  ETpot := GetParFromShell_ETpot( x );
  if ( ETpot < cMin_ETpot ) or ( ETpot > cMax_ETpot ) then begin
    IErr := cInvld_ETpot; Exit;
  end;

  fXG := GetParFromShell_fXG( x );
  if ( fXG < cMin_fXG ) or ( fXG > cMax_fXG ) then begin
    IErr := cInvld_fXG; Exit;
  end;

  cQ := GetParFromShell_cQ( x );
  if ( cQ < cMin_cQ ) or ( cQ > cMax_cQ ) then begin
    IErr := cInvld_cQ; Exit;
  end;

  SL := GetParFromShell_SL( x );
  if ( SL < cMin_SL ) or ( SL > cMax_SL ) then begin
    IErr := cInvld_SL; Exit;
  end;

  Hp := GetParFromShell_Hp( x );
  if ( Hp < cMin_Hp ) or ( Hp > cMax_Hp ) then begin
    IErr := cInvld_Hp; Exit;
  end;

  Wd1 := GetParFromShell_Wd1( x );
  if ( Wd1 < cMin_Wd1 ) or ( Wd1 > cMax_Wd1 ) then begin
    IErr := cInvld_Wd1; Exit;
  end;

  Wd2 := GetParFromShell_Wd2( x );
  if ( Wd2 < cMin_Wd2 ) or ( Wd2 > cMax_Wd2 ) then begin
    IErr := cInvld_Wd2; Exit;
  end;

  Wd3 := GetParFromShell_Wd3( x );
  if ( Wd3 < cMin_Wd3 ) or ( Wd3 > cMax_Wd3 ) then begin
    IErr := cInvld_Wd3; Exit;
  end;

  Wi1 := GetParFromShell_Wi1( x );
  if ( Wi1 < cMin_Wi1 ) or ( Wi1 > cMax_Wi1 ) then begin
    IErr := cInvld_Wi1; Exit;
  end;

  Wi2 := GetParFromShell_Wi2( x );
  if ( Wi2 < cMin_Wi2 ) or ( Wi2 > cMax_Wi2 ) then begin
    IErr := cInvld_Wi2; Exit;
  end;

  Wi3 := GetParFromShell_Wi3( x );
  if ( Wi3 < cMin_Wi3 ) or ( Wi3 > cMax_Wi3 ) then begin
    IErr := cInvld_Wi3; Exit;
  end;

  BD1 := GetParFromShell_BD1( x );
  if ( BD1 < cMin_BD1 ) or ( BD1 > cMax_BD1 ) then begin
    IErr := cInvld_BD1; Exit;
  end;

  BD2 := GetParFromShell_BD2( x );
  if ( BD2 < cMin_BD2 ) or ( BD2 > cMax_BD2 ) then begin
    IErr := cInvld_BD2; Exit;
  end;

  BD3 := GetParFromShell_BD3( x );
  if ( BD3 < cMin_BD3 ) or ( BD3 > cMax_BD3 ) then begin
    {ShowMessage( 'BD3= ' + FloatToSTr( BD3 ) );}
    IErr := cInvld_BD3; Exit;
  end;

  cW := GetParFromShell_cW( x );
  if ( cW < cMin_cW ) or ( cW > cMax_cW ) then begin
    IErr := cInvld_cW; Exit;
  end;

  {###}

//  DefaultPar1 := GetDefaultPar1( KeyValue1 );
//  if ( DefaultPar1 < cMin_ParValue1 ) or ( DefaultPar1 > cMax_ParValue1 ) then begin
//    IErr := cInvld_DefaultPar1; Exit;
//  end;

  {###}

//  DefaultPar2 := GetDefaultPar2( KeyValue1, KeyValue2 );
//  if ( DefaultPar2 < cMinParValue2 ) or ( DefaultPar2 > cMaxParValue2 ) then begin
//    IErr := cInvld_DefaultPar2; Exit;
//  end;

  {###}

//  CalcPar1 := ###
//  if (CalcPar1 < cMinCalcPar) or ###

  {###}

  Result := True; IErr := cNoError;
end; {-Function SetKeyAndParValues}

Function Replace_InitialValues_With_ShellValues( var IErr: Integer): Boolean;
  {-Als de Shell 1-of meer initiele waarden aanlevert voor de array met afhankelijke
    variabelen ('y'), dan kunnen deze waarden hier op deze array worden geplaatst en
    gecontroleerd}
  Function GetParFromShell_PHIT_init( const x: Double ): Double;
  begin
    with EP[ indx-1 ].xDep do
      Result := Items[ cTb_PHIT_init ].EstimateY( 0, Direction );
  end;
begin
  IErr := cNoError; Result := True;
  PHIT_init := GetParFromShell_PHIT_init( x );
  if ( PHIT_init < cMin_PHIT_init ) or ( PHIT_init > cMax_PHIT_init ) then begin
    IErr := cInvld_PHIT_init; Exit;
  end;
  y[ c_PHIT ] := PHIT_init;
  y[ c_dG ] := func_dG( SL, y[ c_PHIT ] );
  y[ c_dV ] := func_dVeq_dG( y[ c_dG ], psi_ae, theta_s, b );
  if ( y[ c_dG ] < psi_ae ) then
    dV_reken := 0
  else
    dV_reken := y[ c_dV ]; {-Trigger wordt nml niet uitgevoerd bij modelinitialisatie }
{$ifdef test} Writeln( f, 'ParFromShell_PHIT_init = ' + FloatToStrF( y[ c_PHIT ], ffFixed, 2, 8 ) );{$endif}
end; {-Replace_InitialValues_With_ShellValues}

Function func_Qdr_j( const PHIT,
                           Hp,
                           BD,
                           Wd,
                           Wi: Double ): Double;
var
  Lj,
  Wj: Double;
  Function Qdr_j( const PHIT, Lj, Wj: Double ): Double;
  begin
    Result := ( PHIT - Lj ) / Wj;
  end;
begin
  Result := 0;
  if ( ( PHIT >= BD ) and ( PHIT >= Hp ) ) then begin {-Drainage if phreatic level higher then bottom of ditch and polder level}
    Lj := max( Hp, BD );
    Wj := Wd;
    Result := Qdr_j( PHIT, Lj, Wj );
  end else if ( ( PHIT <= Hp ) and ( Hp >= BD ) ) then begin {-Infiltration if polder level higher then bottom of ditch and phreatic level smaller then polder level}
    Lj := Hp;
    Wj := Wi;
    Result := Qdr_j( max( PHIT, BD), Lj, Wj );   {-If PHIT<BD, limit the head difference driving the infiltration}
  end else
    Result := 0;
end;

begin {-Procedure DerivsProc}
  {$ifdef test} AssignFile( f, 'testDSmodel204.log' ); Rewrite( f );{$endif}

  for i := 1 to cNrOfDepVar do {-Default speed = 0}
    dydx[ i ] := 0;

  IErr := cUnknownError;

  {-Geef de aanroepende procedure een handvat naar het ModelProfiel}
  if ( nDC > 0 ) then
    aModelProfile := @ModelProfile
  else
    aModelProfile := NIL;

  if not SetKeyAndParValues( IErr ) then
    exit;

  if ( Context = UpdateYstart ) then begin {-Run fase 1}

    {-Optioneel: initiele waarden vervangen door Shell-waarden}
    if not Replace_InitialValues_With_ShellValues( IErr ) then
	  Exit;

    {-Bij Shell-gebruik van het model (indx = cBoot2) dan kan het wenselijk zijn de tijd-as
	  van alle Shell-gegevens te converteren, bijvoorbeeld naar jaren}
//      ### if ( indx = cBoot2 ) then
//        ScaleTimesFromShell( cFromDayToYear, EP ); ###

    IErr := cNoError;
  end else begin {-Run fase 2}   {-Bereken de array met afgeleiden 'dydx'}

   {$ifdef test} Writeln( f, 'PHIT = ' + FloatToStr( y[ c_PHIT ] ) );{$endif}
    with ModelProfile do begin
      y[ c_PHIT ] := DCfunc( Discontinuity, y[ c_PHIT ], GT, Hp,  Context, cDCfunc_Hp, Triggered );
      y[ c_PHIT ] := DCfunc( Discontinuity, y[ c_PHIT ], GT, BD1, Context, cDCfunc_BD1, Triggered );
      y[ c_PHIT ] := DCfunc( Discontinuity, y[ c_PHIT ], GT, BD2, Context, cDCfunc_BD2, Triggered );
      y[ c_PHIT ] := DCfunc( Discontinuity, y[ c_PHIT ], GT, BD3, Context, cDCfunc_BD3, Triggered );
    end;

    {-y[ c_hQ ]: Level quickflow reservoirLevel quickflow reservoir}
    with ModelProfile do begin
      y[ c_hQ ] := DCfunc( DomainBoundary, y[ c_hQ ], LT, 0, Context, cDCfunc_hQ, Triggered );
    end;

    {-Discontinuiteiten in functie W_dV}
    with ModelProfile do begin
      y[ c_dV ] := DCfunc( Discontinuity, y[ c_dV ], Gt, cW, Context, cDCfunc_W_dV_1, Triggered );
    end;
    dV_reken := y[ c_dV ];

    {-Discontinuiteiten in functie dVeq_dG Functie Equilibrium storage deficit (m)
     in afhankelijkheid van de grondwaterstand dG}
    with ModelProfile do begin
      y[ c_dG  ] := DCfunc( Discontinuity, y[ c_dG  ], LT, psi_ae, Context, cDCfunc_dG_psi_ae, Triggered );
      if Triggered then
        dV_reken := 0;
      y[ c_dG  ] := DCfunc( Discontinuity, y[ c_dG  ], LT, 0, Context, cDCfunc_dG_0, Triggered );
      if Triggered then
        dV_reken := 0;
    end;

    dVeq := func_dVeq_dG( y[ c_dG  ], psi_ae, theta_s, b ); {-Equilibrium storage deficit (m)}
    W    := func_W_dV( dV_reken, cW );                      {-Wetness index (-)}
    PV   := P * ( 1 - W );                                  {-Precipitation into vadose zone = P.(1-W) (m/d)}
    PQ   := P * W;                                          {-Precipitation into quickflow reservoir = P.W (m/d)}
    Beta := func_beta_dV( dV_reken, zeta1, zeta2 );         {-Evapotranspiration reduction factor (-)}
    Eact := Beta * ETpot;                                   {-Actual evapotranspiration = Beta.ETpot (m)}

    fQs  := func_fQS( y[ c_hQ ], cQ );                      {-Quickflow (m/d)}
    Q1 := func_Qdr_j( y[ c_PHIT ], Hp, BD1, Wd1, Wi1 );     {-Discharge of groundwater to/from primary drainage system (m/d)}
    Q2 := func_Qdr_j( y[ c_PHIT ], Hp, BD2, Wd2, Wi2 );     {-Discharge of groundwater to/from secundary drainage system (m/d)}
    Q3 := func_Qdr_j( y[ c_PHIT ], Hp, BD3, Wd3, Wi3 );     {-Discharge of groundwater to/from tertiary drainage system (m/d)}
    fGS := Q1+Q2+Q3;                                        {-Discharge of groundwater to/from drainage system (=Q1+Q2+Q3) (m/d)}

    {-Change of state variables}
    dydx[ c_dV ]   := -( fXG + PV - Eact - fGS );           {-Storage deficit (m)}
    if ( dV_reken > 0 ) then
      dydx[ c_dG ] := ( y[ c_dV ] - dVeq ) / cV
    else
      dydx[c_dG ] := dydx[ c_dV ];                          {-Groundwater depth (m-mv)}

    dydx[ c_PHIT ] := -dydx[ c_dG ];

{$ifdef test} Writeln( f, 'PHIT = ' + FloatToStrF( y[ c_PHIT ], ffFixed, 5, 8 ) );{$endif}
{$ifdef test} Writeln( f, 'dydx[ c_PHIT ]   = ' + FloatToStrF( dydx[ c_PHIT ], ffFixed, 5, 8 ) );{$endif}
{$ifdef test} Writeln( f, 'dV               = ' + FloatToStrF( y[ c_dV ], ffFixed, 2, 8 ) );{$endif}
{$ifdef test} Writeln( f, 'dydx[ c_dV ]     = ' + FloatToStrF( dydx[ c_dV ], ffFixed, 5, 8 ) );{$endif}
{$ifdef test} Writeln( f, 'Q1   = ' + FloatToStrF( Q1, ffFixed, 5, 8 ) );{$endif}
{$ifdef test} Writeln( f, 'Q2   = ' + FloatToStrF( Q2, ffFixed, 5, 8 ) );{$endif}
{$ifdef test} Writeln( f, 'Q3   = ' + FloatToStrF( Q3, ffFixed, 5, 8 ) );{$endif}

    dydx[ c_hQ ]   := PQ - fQs;                             {-hQ: Level quickflow reservoir = P.W (m)}

    {-Change of dependent variables}
    dydx[ c_W ]    := W;       {-W: Wetness index (-)}
    dydx[ c_BdV ]  := Beta;    {-BdV: Evapotranspiration reduction factor}
    dydx[ c_dVeq ] := dVeq;    {-dVeq: Equilibrium storage deficit (m)}

    {-External fluxes: input}
    dydx[ c_P ]      := P;     {-(=PQ+PV): Precipitation (m/d)}
    dydx[ c_ETpot  ] := ETpot; {-Potential evapotranspiration (m/d)}
    dydx[ c_fXG  ]   := fXG;   {-Seepage (up/down)/extraction (m/d)}

    {-External fluxes: output}
    dydx[ c_Eact ] := Eact;    {-Eact: Actual evapotranspiration (m/d)}
    dydx[ c_Q1 ]   := Q1;      {-Q1: Discharge of groundwater to/from primary drainage system (m/d)}
    dydx[ c_Q2 ]   := Q2;      {-Q2: Discharge of groundwater to/from secundary drainage system (m/d)}
    dydx[ c_Q3 ]   := Q3;      {-Q3: Discharge of groundwater to/from tertiary drainage system (m/d)}
    dydx[ c_fGS ]  := fGS;     {-fGS: Discharge of groundwater to/from drainage system (=Q1+Q2+Q3) (m/d)}
    dydx[ c_fQS ]  := fQS;     {-fQS: Quickflow (m/d)}
    {-Internal fluxes}
    dydx[ c_PV ]   := PV;      {-PV: Precipitation into vadose zone = P.(1-W) (m/d)}
    dydx[ c_PQ ]   := PQ;      {-PQ: Precipitation into quickflow reservoir = P.W (m/d)}
    {-Other}
    dydx[ c_Err ] := ( fXG + P ) - ( Eact + fGS + fQS ) - ( -dydx[ c_dV ] + dydx[ c_hQ ]);

  end;

  {$ifdef test} CloseFile( f );{$endif}
end; {-DerivsProc}

Function DefaultBootEP( const EpDir: String; const BootEpArrayOption: TBootEpArrayOption; var EP: TExtParArray ): Integer;
  {-Initialiseer de meest elementaire gegevens van het model. Shell-gegevens worden door deze
    procedure NIET verwerkt}
Procedure SetMinMaxKeyAndParValues;
begin
  with EP[ cEP0 ].xInDep.Items[ cTb_MinMaxValKeys ] do begin
    cMin_P     :=        GetValue( 1, 1 );  {rij, kolom}
    cMax_P     :=        GetValue( 1, 2 );
    cMin_ETpot :=        GetValue( 2, 1 );
    cMax_ETpot :=        GetValue( 2, 2 );
    cMin_fXG   :=        GetValue( 3, 1 );
    cMax_fXG   :=        GetValue( 3, 2 );
    cMin_Soil  := Trunc( GetValue( 4, 1 ) );
    cMax_Soil  := Trunc( GetValue( 4, 2 ) );
    cMin_cQ    :=        GetValue( 5, 1 );
    cMax_cQ    :=        GetValue( 5, 2 );
    cMin_SL    :=        GetValue( 6, 1 );
    cMax_SL    :=        GetValue( 6, 2 );
    cMin_Hp    :=        GetValue( 7, 1 );
    cMax_Hp    :=        GetValue( 7, 2 );
    cMin_Wd1   :=        GetValue( 8, 1 );
    cMax_Wd1   :=        GetValue( 8, 2 );
    cMin_Wd2   :=        GetValue( 9, 1 );
    cMax_Wd2   :=        GetValue( 9, 2 );
    cMin_Wd3   :=        GetValue( 10, 1 );
    cMax_Wd3   :=        GetValue( 10, 2 );
    cMin_Wi1   :=        GetValue( 11, 1 );
    cMax_Wi1   :=        GetValue( 11, 2 );
    cMin_Wi2   :=        GetValue( 12, 1 );
    cMax_Wi2   :=        GetValue( 12, 2 );
    cMin_Wi3   :=        GetValue( 13, 1 );
    cMax_Wi3   :=        GetValue( 13, 2 );
    cMin_BD1   :=        GetValue( 14, 1 );
    cMax_BD1   :=        GetValue( 14, 2 );
    cMin_BD2   :=        GetValue( 15, 1 );
    cMax_BD2   :=        GetValue( 15, 2 );
    cMin_BD3   :=        GetValue( 16, 1 );
    cMax_BD3   :=        GetValue( 16, 2 );
    cMin_PHIT_init :=    GetValue( 17, 1 );
    cMax_PHIT_init :=    GetValue( 17, 2 );
    cMin_cW :=           GetValue( 18, 1 );
    cMax_cW :=           GetValue( 18, 2 );
  end;
end;
Begin
  Result := DefaultBootEPFromTextFile( EpDir, BootEpArrayOption, cModelID, cNrOfDepVar, nDC, cNrXIndepTblsInEP0,
                                       cNrXdepTblsInEP0, Indx, EP );
  if ( Result = cNoError ) then begin
    SetMinMaxKeyAndParValues;
    {###SetAnalytic_DerivsProc( True, EP );} {-Ref. 'USpeedProc.pas'}
  end;
end;

Function TestBootEP( const EpDir: String; const BootEpArrayOption: TBootEpArrayOption; var EP: TExtParArray ): Integer;
  {-Deze boot-procedure verwerkt alle basisgegevens van het model en leest de Shell-gegevens
    uit een bestand. Na initialisatie met deze boot-procedure is het model dus gereed om
	'te draaien'. Deze procedure kan dus worden gebruikt om het model 'los' van de Shell te
	testen}
Begin
  Result := DefaultBootEP( EpDir, BootEpArrayOption, EP );
  if ( Result <> cNoError ) then
    exit;
  Result := DefaultTestBootEPFromTextFile( EpDir, BootEpArrayOption, cModelID, cnRP + cnSQ + cnRQ, Indx, EP );
  if ( Result <> cNoError ) then
    exit;
  SetReadyToRun( EP);
end;

Function BootEPForShell( const EpDir: String; const BootEpArrayOption: TBootEpArrayOption; var EP: TExtParArray ): Integer;
  {-Deze procedure maakt het model gereed voor Shell-gebruik.
    De xDep-tables in EP[ indx-1 ] worden door deze procedure NIET geinitialiseerd omdat deze
	gegevens door de Shell worden verschaft }
begin
  Result := DefaultBootEP( EpDir, cBootEPFromTextFile, EP );
  if ( Result = cNoError ) then
    Result := DefaultBootEPForShell( cnRP, cnSQ, cnRQ, Indx, EP );
end;

Exports DerivsProc       index cModelIndxForTDSmodels, {999}
        DefaultBootEP    index cBoot0, {1}
        TestBootEP       index cBoot1, {2}
        BootEPForShell   index cBoot2; {3}

begin
  {-Dit zgn. 'DLL-Main-block' wordt uitgevoerd als de DLL voor het eerst in het geheugen wordt
    gezet (Reason = DLL_PROCESS_ATTACH)}
  DLLProc := @MyDllProc;
  Indx := cBootEPArrayVariantIndexUnknown;
  if ( nDC > 0 ) then
    ModelProfile := TModelProfile.Create( nDC );
end.
