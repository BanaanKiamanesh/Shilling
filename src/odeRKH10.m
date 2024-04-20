function [Time, Y] = odeRKH10(f, TSpan, Y0, h)
    %ODERKH10 10th Order Runge-Kutta Hairer's Method Implementation
    %
    % Method Properties:
    %     Method Name:
    %                         10th order Runge-Kutta Hairer
    %     Introduced in Year:
    %                         1978
    %     Method Type:
    %                         Fixed Time Step, Explicit
    %     Order:
    %                         10
    %     Number of Stages:
    %                         17
    %     Number of Registers:
    %                         17
    %     Links:
    %         1. https://www.researchgate.net/publication/31221486_A_Runge-Kutta_Method_of_Order_10
    %
    %
    % Inputs:
    %   f: function handle for the ODE
    %   TSpan: time span as [t0, tf]
    %   Y0: initial condition as a column vector
    %   h: step size (default: 0.01)
    %
    %
    % Outputs:
    %   Time: time vector associated with the integration
    %   Y: solved ode state evolution matrix
    %
    %
    % Example Usage:
    %   f = @(t, x) [x(2); (1 - x(1)^2)*x(2) - x(1)]; % Van der Pol ODE
    %   TSpan= [0, 20];    % Time span
    %   Y0= [2; 0];        % Initial condition
    %   h = 0.01;          % Initial step size
    %   [T, Y] = odeRKH10(f, TSpan, Y0, h);
    %
    %   % Plot results
    %   figure;
    %   plot(T, Y, 'LineWidth', 2);
    %   xlabel('Time [s]');
    %   ylabel('State Variables');
    %   title('Van der Pol Equation Solution');
    %   grid on;
    %   axis([0, 20, -3, 3]);
    %
    % Reference
    %   * Ernst Hairer, "A Runge-Kutta Method of Order 10",
    %       January 1978, IMA Journal of Applied Mathematics 21(1)

    % Set default values if not provided
    if nargin < 4
        h = 0.01;
    end

    t0 = TSpan(1);
    tf = TSpan(2);

    % Initial setup
    t = t0;
    y = Y0;

    % Method Params
    digits(85);
    a2  = vpa(0.5233584004620047139632937023215170497515953383496781610502942213792083195573343614542);
    a3  = vpa(0.5265091001416125727329516734775408259523008085442588621240322448552569517961160776999);
    a4  = vpa(0.7897636502124188590994275102163112389284512128163882931860483672828854276941741165498);
    a5  = vpa(0.3939235701256720143227738119996521367488350094732736567444747389961242969755195414769);
    a6  = vpa(0.7666539862535505911932668693560686601417881683220066628197494388337786524595845606945);
    a7  = vpa(0.2897636502124188590994275102163112389284512128163882931860483672828854276941741165498);
    a8  = vpa(0.1084776892195672933536461100396272126536495082872530661076180009074872219675135232227);
    a9  = vpa(0.3573842417596774518429245029795604640404982636367873040901247917361510345429002009092);
    a10 = vpa(0.8825276619647323464255014869796690751828678442680521196637911779185276585194132570617);
    a11 = vpa(0.6426157582403225481570754970204395359595017363632126959098752082638489654570997990908);
    a12 = vpa(0.1174723380352676535744985130203309248171321557319478803362088220814723414805867429383);
    a13 = vpa(0.7666539862535505911932668693560686601417881683220066628197494388337786524595845606945);
    a14 = vpa(0.2897636502124188590994275102163112389284512128163882931860483672828854276941741165498);
    a15 = vpa(0.5265091001416125727329516734775408259523008085442588621240322448552569517961160776999);
    a16 = vpa(0.5233584004620047139632937023215170497515953383496781610502942213792083195573343614542);

    b21   = vpa(0.5233584004620047139632937023215170497515953383496781610502942213792083195573343614542);
    b31   = vpa(0.2616697163778127283312402097548997641973627614039327706102102491953717704187699219879);
    b32   = vpa(0.2648393837637998444017114637226410617549380471403260915138219956598851813773461557120);
    b41   = vpa(0.1974409125531047147748568775540778097321128032040970732965120918207213569235435291375);
    b43   = vpa(0.5923227376593141443245706326622334291963384096122912198895362754621640707706305874124);
    b51   = vpa(0.1973205486287023067036649485978952117578825584337199991656132317825743833021058860612);
    b53   = vpa(0.2950833340926721918228255598274359230145509784596048092593866299668586683096812321792);
    b54   = vpa(-0.9848031259570248420371669642567899802359852742005115168052512275330875463626757676351e-1);
    b61   = vpa(0.1313134173444616536130177999345909470542768366990469474228776563495603985387834598888);
    b64   = vpa(0.1101544395386396206773696967716892932905881833462590372574439152868807925960672597269);
    b65   = vpa(0.5251861293704493169028793726497884197969231482767006781394278671973374613247338410788);
    b71   = vpa(0.1342003418463226002727476951680931444018781918996634475042938727232827990267510411664);
    b74   = vpa(0.6960887032881160802299824047678314828416281106463280160258778625630871045892118300249);
    b75   = vpa(0.2504977215703398097125518092509800218006823479445679226772611578145189247821334145350);
    b76   = vpa(-0.7910231164923596311158543989705934101157374376741710930213845258180034007039221691766);
    b81   = vpa(0.7221827418966261942005081845517049770474117717860435538167259617193885317787239789517e-1);
    b85   = vpa(-0.5833632293645610716380606138930651874161115931850840194371530896953184153028603826158e-1);
    b86   = vpa(0.3047557668574525220174950070294036092519082552287015783049982467857131629284936232933e-2);
    b87   = vpa(0.9154818029778625587722640290346919759800040787487009688661073123722307869064222735619e-1);
    b91   = vpa(0.3125500813516617947050120528263476612150928811800844295622424453129979711857118942140e-1);
    b96   = vpa(0.1091238215424128929294834955207716494619546474783251185719596191189550956073995218558e-3);
    b97   = vpa(0.1567257586309938356246107465648794708917441337700138495832023584660674808897051732313);
    b98   = vpa(0.1692943511719750238548830676365254553777828871012866864321262291196648014390164387346);
    b101  = vpa(0.1190660441466861924216884258080925099509546061156163781761478570338134316043617322173e-1);
    b106  = vpa(0.2834370820246027860255992266982188665428702252125274101591962479489267620092644600260);
    b107  = vpa(-0.4163121675706282353724276181356300212164192944619403444080980777942870933794472685216);
    b108  = vpa(0.2646463339497663668210902091085361027053564926326924812994390366777834273576276339310);
    b109  = vpa(0.7388498091463228097090708267277348761559649602732109347956391853827232193715322584046);
    b111  = vpa(0.2340657369133197891470838377984007842503946857756845416362339865999662377057916960290e-1);
    b116  = vpa(0.9449313018949365401300253095605614324982516623777346096448800625130954018295939047740e-1);
    b117  = vpa(-0.2728720559019952606363092580665963250433705067252372208829562550632597611313198545757);
    b118  = vpa(0.2240220461156057997944315522518131846261124699330640294272458923970695162204120486767);
    b119  = vpa(0.6043814410751657569719347222576085340011863610739072982875214128849988434755301302413);
    b1110 = vpa(-0.3081537692927938090069243415828207929929122273386332605004724686626579706106108533173e-1);
    b121  = vpa(0.4544377531017616315765389908153096498645890941991961177836633137902353666760633687088e-1);
    b126  = vpa(-0.1187996671864028586765254219285356343376285990176386474891150929245660355742130860183e-2);
    b127  = vpa(0.1203565499092261097966188217234362058515446695047694116231895919296441480090125575596e-1);
    b128  = vpa(0.7512690298764966821627521371565572140275315500656240513504528112598150181393445048666e-1);
    b129  = vpa(-0.1822092409888012403141186105974838892761579859192182074696828369088959767419077567178e-1);
    b1210 = vpa(-0.2571528540841043468806376221771396205460354181513400383106090430345956162981031278207e-3);
    b1211 = vpa(0.4532078371347468185965270952011502734303744355238469520648294046672741844375709484540e-2);
    b131  = vpa(0.1767137782592772030958798765711993346076326211800572275450227165783753236705910865492);
    b134  = vpa(0.1101544395386396206773696967716892932905881833462590372574439152868807925960672597269);
    b135  = vpa(0.5251861293704493169028793726497884197969231482767006781394278671973374613247338410788);
    b136  = vpa(-0.4716207672801957948798217912152359376250630852495511063738116933651587031904328351457);
    b137  = vpa(0.8990310498491875266368990071875152922763468480002185650326986125011485318362907529907);
    b138  = vpa(-0.7467230306916289638599602008088168117750310724922743198498253813592425510843163068237);
    b139  = vpa(-10.017101516756146040853186972006065972987027196800421553809421717321497529906933631477);
    b1310 = vpa(0.1263508715195988962951307827687648346421985369266969430473204298972536422365713122404);
    b1311 = vpa(0.5660138272355064270682732249907470012763799581315503842554078250210353407723389384909);
    b1312 = vpa(0.5986492052088624001098038724464832066388402270027708075754868643976463442046741430643);
    b141  = vpa(0.1277534947480869822694777006880571541639616513225826576695303067404023367054772185702);
    b144  = vpa(0.6960887032881160802299824047678314828416281106463280160258778625630871045892118300249);
    b145  = vpa(0.2504977215703398097125518092509800218006823479445679226772611578145189247821334145350);
    b146  = vpa(-0.7368246436028416867609246757454535374296880219263938462439002090823944915566264811824);
    b147  = vpa(-0.2778578777108241826773273374900723250222301109862216853553157201018147214465526588169);
    b148  = vpa(-0.5997526313598403501296884799197753021563938240370770948150479630779446286262003432092);
    b149  = vpa(0.2024692338910704693500237585621903123505161701229471467587157451308903694383321235511);
    b1410 = vpa(0.5432036982363849780600684652634443601468189969678666775046718813224989883416104871445e-2);
    b1411 = vpa(-0.1074472474155047920101206919894381337125444664272205024314798936418733258920769563370e-1);
    b1412 = vpa(0.6951688484570234004700591858164146072357628221597117426434839740273190245052113679250);
    b1413 = vpa(-0.6246651130952503394431547116755180508600167575701318270645551618021614799102076408562e-1);
    b151  = vpa(0.2616697163778127283312402097548997641973627614039327706102102491953717704187699219879);
    b152  = vpa(0.2648393837637998444017114637226410617549380471403260915138219956598851813773461557120);
    b156  = vpa(-0.1998011270205324791079663580830885049848273745422651189682301346802905866051733476638);
    b157  = vpa(-0.6510499873052827124921914489683813643155863882516440645794556633240216912803403931627);
    b1513 = vpa(0.1998011270205324791079663580830885049848273745422651189682301346802905866051733476638);
    b1514 = vpa(0.6510499873052827124921914489683813643155863882516440645794556633240216912803403931627);
    b161  = vpa(0.5233584004620047139632937023215170497515953383496781610502942213792083195573343614542);
    b163  = vpa(-0.5558812136754302060726143105309293455559184141943321053532734480099926250948077261183);
    b1615 = vpa(0.5558812136754302060726143105309293455559184141943321053532734480099926250948077261183);
    b171  = vpa(0.5732079543206559103114261705103983656495216504867462310285994428078568043160654439795e-1);
    b172  = vpa(-0.5499710763899945608115841896290187887481592249811405834035066676393750158953834290913);
    b173  = vpa(-0.6499374174008749135116607420010890619711618624173024222960650740195874521599402439688);
    b176  = vpa(-10.061667370401756207240019539023157074172524666307437022389776456477183230723296269940);
    b177  = vpa(-0.4040156689806358294269682234212183308262562023912486365220642577870402491555711062480e-1);
    b178  = vpa(-0.1828302366407607254710272774065261039379052622607190097473388370699414811305446343873);
    b179  = vpa(-0.3336592706492786845666575661828162687906558601961826440714525336287466822150370633233);
    b1710 = vpa(0.3956485423760567568801345107166015519577734440834727480004748180136901286634710478955);
    b1711 = vpa(0.6950570494599735891002099282005158129027126868215679095299345058137097320818106877162);
    b1712 = vpa(0.2714873764573748588377263058539220945263829691804714618529052530298982146739754552950);
    b1713 = vpa(0.6071810560414041202873774349794680164722661545496003750296400378855628528787164400954);
    b1714 = vpa(0.5918636248229842840838104081530739675596239893196764223449596939309288102548549028752);
    b1715 = vpa(0.6499374174008749135116607420010890619711618624173024222960650740195874521599402439688);
    b1716 = vpa(0.5499710763899945608115841896290187887481592249811405834035066676393750158953834290913);

    c1  = vpa(0.3333333333333333333333333333333333333333333333333333333333333333333333333333333333333e-1);
    c2  = vpa(-0.3846153846153846153846153846153846153846153846153846153846153846153846153846153846154e-1);
    c3  = vpa(-0.9090909090909090909090909090909090909090909090909090909090909090909090909090909090909e-1);
    c6  = vpa(-0.1348314606741573033707865168539325842696629213483146067415730337078651685393258426966);
    c7  = vpa(-0.1111111111111111111111111111111111111111111111111111111111111111111111111111111111111);
    c9  = vpa(0.2774291885177431765083602625606543404285043197180408363394722409866844803871713937960);
    c10 = vpa(0.1892374781489234901583064041060123262381623469486258303271944256799821862794952728707);
    c11 = vpa(0.2774291885177431765083602625606543404285043197180408363394722409866844803871713937960);
    c12 = vpa(0.1892374781489234901583064041060123262381623469486258303271944256799821862794952728707);
    c13 = vpa(0.1348314606741573033707865168539325842696629213483146067415730337078651685393258426966);
    c14 = vpa(0.1111111111111111111111111111111111111111111111111111111111111111111111111111111111111);
    c15 = vpa(0.9090909090909090909090909090909090909090909090909090909090909090909090909090909090909e-1);
    c16 = vpa(0.3846153846153846153846153846153846153846153846153846153846153846153846153846153846154e-1);
    c17 = vpa(0.3333333333333333333333333333333333333333333333333333333333333333333333333333333333333e-1);

    % Preallocate arrays to store values
    StepNum = ceil((tf - t0)  /  h);
    Time    = zeros(StepNum, 1);
    Y       = zeros(StepNum, length(Y0));
    Time(1) = t;
    Y(1, :) = y';

    % Main Loop
    Idx = 1;
    while t < tf
        % Method
        k1  = f(        t,                                                                                                                                                               y);
        k2  = f( t + a2*h,                                                                                                                                                  y + h*(b21*k1));
        k3  = f( t + a3*h,                                                                                                                                         y + h*(b31*k1 + b32*k2));
        k4  = f( t + a4*h,                                                                                                                                         y + h*(b41*k1 + b43*k3));
        k5  = f( t + a5*h,                                                                                                                                y + h*(b51*k1 + b53*k3 + b54*k4));
        k6  = f( t + a6*h,                                                                                                                                y + h*(b61*k1 + b64*k4 + b65*k5));
        k7  = f( t + a7*h,                                                                                                                       y + h*(b71*k1 + b74*k4 + b75*k5 + b76*k6));
        k8  = f( t + a8*h,                                                                                                                       y + h*(b81*k1 + b85*k5 + b86*k6 + b87*k7));
        k9  = f( t + a9*h,                                                                                                                       y + h*(b91*k1 + b96*k6 + b97*k7 + b98*k8));
        k10 = f(t + a10*h,                                                                                                         y + h*(b101*k1 + b106*k6 + b107*k7 + b108*k8 + b109*k9));
        k11 = f(t + a11*h,                                                                                             y + h*(b111*k1 + b116*k6 + b117*k7 + b118*k8 + b119*k9 + b1110*k10));
        k12 = f(t + a12*h,                                                                                 y + h*(b121*k1 + b126*k6 + b127*k7 + b128*k8 + b129*k9 + b1210*k10 + b1211*k11));
        k13 = f(t + a13*h,                                                 y + h*(b131*k1 + b134*k4 + b135*k5 + b136*k6 + b137*k7 + b138*k8 + b139*k9 + b1310*k10 + b1311*k11 + b1312*k12));
        k14 = f(t + a14*h,                                     y + h*(b141*k1 + b144*k4 + b145*k5 + b146*k6 + b147*k7 + b148*k8 + b149*k9 + b1410*k10 + b1411*k11 + b1412*k12 + b1413*k13));
        k15 = f(t + a15*h,                                                                                           y + h*(b151*k1 + b152*k2 + b156*k6 + b157*k7 + b1513*k13 + b1514*k14));
        k16 = f(t + a16*h,                                                                                                                           y + h*(b161*k1 + b163*k3 + b1615*k15));
        k17 = f(    t + h, y + h*(b171*k1 + b172*k2 + b173*k3 + b176*k6 + b177*k7 + b178*k8 + b179*k9 + b1710*k10 + b1711*k11 + b1712*k12 + b1713*k13 + b1714*k14 + b1715*k15 + b1716*k16));

        % Update Values
        y = y + h*(c1*k1 + c2*k2 + c3*k3 + c6*k6 + c7*k7 + c9*k9 + c10*k10 + c11*k11 + c12*k12 + c13*k13 + c14*k14 + c15*k15 + c16*k16 + c17*k17);
        t = t + h;

        % Store Values
        Idx       = Idx + 1;
        Time(Idx) = t;
        Y(Idx, :) = y';
    end
end
