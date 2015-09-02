(* ::Package:: *)

(************************************************************************)
(* This file was generated automatically by the Mathematica front end.  *)
(* It contains Initialization cells from a Notebook file, which         *)
(* typically will have the same name as this file except ending in      *)
(* ".nb" instead of ".m".                                               *)
(*                                                                      *)
(* This file is intended to be loaded into the Mathematica kernel using *)
(* the package loading commands Get or Needs.  Doing so is equivalent   *)
(* to using the Evaluate Initialization Cells menu command in the front *)
(* end.                                                                 *)
(*                                                                      *)
(* DO NOT EDIT THIS FILE.  This entire file is regenerated              *)
(* automatically each time the parent Notebook file is saved in the     *)
(* Mathematica front end.  Any changes you make to this file will be    *)
(* overwritten.                                                         *)
(************************************************************************)



ClearAll[cn,\[CapitalDelta],l,n,\[Epsilon],coeffs];






SetNDer[nzder_,prec_,zv_]:=
(* this function computs coeff[i,j] with i,j<=nzder *)
Module[{\[Rho]0,\[Rho]der,F,coeff,coeffTab,zdermax,f,zz,z},(* \[Rho]der[n] is the n-th derivative of rho(z) at z=1/2 *)
$MaxExtraPrecision=100;
Clear[\[Rho]0,\[Rho]der];
t1=AbsoluteTiming[\[Rho]0[0]=(2-z-2Sqrt[1-z])/z;
\[Rho]0[n_]:=\[Rho]0[n]=D[\[Rho]0[n-1],z];
\[Rho]der[n_]:=\[Rho]der[n]=((\[Rho]0[n]/.z->zv//Simplify))//Expand//Simplify][[1]];

ClearAll[F];

coeff[i_,j_]:=BellY[i,j,Table[\[Rho]der[p],{p,1,i-j+1}]];
coeffTab = Table[coeff[i,j],{i,nzder},{j,i}];
zdermax = nzder;
Return[SetPrecision[coeffTab,prec]];
];





Options[Recurse]={coeffTab->None};
Recurse[eps_,ll_,prec_,mmax_,nmax_,kmax_,OptionsPattern[]]:=Module[{\[Epsilon],l,fudge,expr,zdermax,t1,t2,t3,t4,coeffs,cn,cnRhoDer,CRec,cs,\[Rho]num,coeffTb,z,casCoeffs},(

ClearAll[cn,cnRhoDer,CRec];

\[Epsilon]=eps;l=ll;
fudge=10^-50;

coeffs[n_]:={((-3+n+\[CapitalDelta]) (-5-l+2 n+\[CapitalDelta]) (-2+n+\[Epsilon]) (-5+l+2 n+\[CapitalDelta]+2 \[Epsilon]))/(n (-1+l+2 n+\[CapitalDelta]) (-1-l+2 n+\[CapitalDelta]-2 \[Epsilon]) (-1+n+\[CapitalDelta]-\[Epsilon])),1/(4 n (-1+l+2 n+\[CapitalDelta]) (-1-l+2 n+\[CapitalDelta]-2 \[Epsilon]) (-1+n+\[CapitalDelta]-\[Epsilon])) (-3 l (-1+\[CapitalDelta]) (-1+\[CapitalDelta]-2 \[Epsilon]) (l+2 \[Epsilon])-(-4+2 n+\[CapitalDelta]) (-166+24 n^3+3 \[CapitalDelta]^3+4 n^2 (9 \[CapitalDelta]+4 (-8+\[Epsilon]))+4 \[CapitalDelta]^2 (-8+\[Epsilon])+78 \[Epsilon]+4 \[Epsilon]^2+\[CapitalDelta] (123-34 \[Epsilon]-4 \[Epsilon]^2)+2 n (123+9 \[CapitalDelta]^2+8 \[CapitalDelta] (-8+\[Epsilon])-34 \[Epsilon]-4 \[Epsilon]^2))+(3+12 (-2+n)^2-2 \[Epsilon]+4 (-2+n) (1+3 \[CapitalDelta]+\[Epsilon])+\[CapitalDelta] (2+3 \[CapitalDelta]+2 \[Epsilon])) (l^2+2 l \[Epsilon]+\[CapitalDelta] (\[CapitalDelta]-2 (1+\[Epsilon])))),((-2+2 n+\[CapitalDelta]) ((-2+2 n+\[CapitalDelta]) (31+12 n^2-16 \[CapitalDelta]+3 \[CapitalDelta]^2+4 n (-8+3 \[CapitalDelta]))-2 (9+8 n^2-7 \[CapitalDelta]+2 \[CapitalDelta]^2+2 n (-7+4 \[CapitalDelta])) \[Epsilon]-4 (-5+2 n+\[CapitalDelta]) \[Epsilon]^2)+3 (-2 (-2+2 n+\[CapitalDelta])+l (-1+\[CapitalDelta]) (-1+\[CapitalDelta]-2 \[Epsilon]) (l+2 \[Epsilon]))+(-19-12 n^2-3 \[CapitalDelta]^2-2 \[Epsilon]+2 \[CapitalDelta] (7+\[Epsilon])+4 n (7-3 \[CapitalDelta]+\[Epsilon])) (l^2+2 l \[Epsilon]+\[CapitalDelta] (\[CapitalDelta]-2 (1+\[Epsilon]))))/(4 n (-1+l+2 n+\[CapitalDelta]) (-1-l+2 n+\[CapitalDelta]-2 \[Epsilon]) (-1+n+\[CapitalDelta]-\[Epsilon]))};

zdermax=2nmax+mmax;


\[Rho]num=SetPrecision[3-2Sqrt[2],prec];

z=1/2; (*If we want to generate blocks at some other point this needs to be changed*)
coeffTb=OptionValue[coeffTab];

casCoeffs[m_,n_,\[Epsilon]_,l_]={-((m (2-3 m+m^2) )/(4 z^2 (-2+2 z))),-(((-1+m) m (-2+6 z) )/(4 z^2 (-2+2 z))),-(m (22+m^2+12 n^2+2 \[Epsilon]-2 n (17+2 \[Epsilon])+m (-13+12 n+2 \[Epsilon])))/(8 z^2 (-2+2 z) (-1+2 n+2 \[Epsilon])),(m (4-6 z) )/(2 z (-2+2 z)),1/(8 z^2 (-2+2 z) (-1+2 n+2 \[Epsilon])) (-2 z (3 m^2+2 (-1+n) (-5+6 n-2 \[Epsilon])+m (-27+24 n+2 (2+2 \[Epsilon])))+2 (m^2+m (-9+8 n)+2 (3+2 n^2-1/2 \[CapitalDelta] (-2+\[CapitalDelta]-2 \[Epsilon])+2 \[Epsilon]-n (5+2 \[Epsilon])-1/2 l (l+2 \[Epsilon])))),-(((-1+n) (-6+3 m+4 n-2 \[Epsilon]) )/(8 z^2 (-2+2 z) (-1+2 n+2 \[Epsilon]))),+((4 (-4+m+4 n)-2 z (-10+3 m+12 n+2 \[Epsilon])))/(4 z (-2+2 z) (-1+2 n+2 \[Epsilon])),-(((-1+n) (-2+6 z))/(8 z^2 (-2+2 z) (-1+2 n+2 \[Epsilon]))),+(1/(6-4 n-2 (2+2 \[Epsilon])))}//Simplify;

casInds[m_,n_]={{-3+m,n},{-2+m,n},{-1+m,-1+n},{-1+m,n}, {m,-1+n},{1+m,-2+n},{1+m,-1+n}, {2+m,-2+n},{2+m,-1+n}};

(*************************)





cn[n_]:=cn[n]=Chop[Apart[Expand[Sum[coeffs[n][[i]]cn[n-4+i],{i,3}]],\[CapitalDelta]]];
cn[0]=1;
cn[x_/;x<0]=0;

(* Recall, for identical scalar case, series in powers of rho^2n, hence the factor of 2 below*)
(* Power series are rho^\Delta Sum( cnrhoder[n,i] rho^(2n-i));  *)
(* fudge factor is required to prevent certain cancellations from occurring; otherwise different components will have different poles *)



cnRhoDer[n_,i_]:=cnRhoDer[n,i]=Expand[(\[CapitalDelta]+2n-(i-1)+fudge)cnRhoDer[n,i-1]]/.{\[CapitalDelta] r[1,a_]:>  1+ a r[1,a],\[CapitalDelta]^2 r[1,a_]:>  \[CapitalDelta]+a+ a^2 r[1,a]}/.{\[CapitalDelta] r[2,a_]:>  a r[2,a]+ r[1,a],\[CapitalDelta]^2 r[2,a_]:>  1+a^2 r[2,a]+2a r[1,a],\[CapitalDelta]^3 r[2,a_]->2a+\[CapitalDelta]+a^3 r[2,a]+3a^2 r[1,a]};
cnRhoDer[n_,0]:=cnRhoDer[n,0]=cn[n]/.{1/\[CapitalDelta]:>r[1,0],1/\[CapitalDelta]^2:>r[2,0],1/(a_+b_ \[CapitalDelta]):> 1/b r[1,-a/b],1/(a_+\[CapitalDelta]):> r[1,-a]}/.{Power[a_+b_ \[CapitalDelta],m_]:> 1/b^-m r[-m,-a/b],Power[a_+\[CapitalDelta],m_]:> r[-m,-a]};

t1=AbsoluteTiming[
Do[cn[k],{k,0,kmax}]];


(*Need to convert to z=zb derivatives (or more precisely the a derivatives *)
t2=AbsoluteTiming[
AllDer[ll,0,0]=(Sum[Expand[cn[k]\[Rho]num^(2k)],{k,0,kmax}]//Expand)/.{1/\[CapitalDelta]:>r[1,0],1/\[CapitalDelta]^2:>r[2,0],1/(a_+b_ \[CapitalDelta]):> 1/b r[1,-a/b],1/(a_+\[CapitalDelta]):> r[1,-a]}/.{Power[a_+b_ \[CapitalDelta],m_]:> 1/b^-m r[-m,-a/b],Power[a_+\[CapitalDelta],m_]:> r[-m,-a]}//Simplify;
];

t3=AbsoluteTiming[
Do[
AllDer[ll,m,0]=(1/2^m Sum[Expand[coeffTb[[m,i]]cnRhoDer[k,i]\[Rho]num^(2k-i)]/.{\[CapitalDelta] r[1,a_]:>  1+ a r[1,a],\[CapitalDelta]^2 r[1,a_]:>  \[CapitalDelta]+a+ a^2 r[1,a]}/.{\[CapitalDelta] r[2,a_]:>  a r[2,a]+ r[1,a],\[CapitalDelta]^2 r[2,a_]:>  1+a^2 r[2,a]+2a r[1,a],\[CapitalDelta]^3 r[2,a_]->2a+\[CapitalDelta]+a^3 r[2,a]+3a^2 r[1,a]},{k,0,kmax},{i,1,Length[coeffTb[[m]]]}]);
,{m,1,zdermax}];
];

(* Use Casimir rec *)
CRec[m_/;m<0,n_]:=0;
CRec[m_,n_/;n<0]:=0;
Do[CRec[m,0]=AllDer[ll,m,0],{m,0,zdermax}];

t4=AbsoluteTiming[
Do[
cs=casCoeffs[m,n,eps,ll]//Simplify;

AllDer[ll,m,n]=CRec[m,n]=Expand[Sum[cs[[i]]CRec@@casInds[m,n][[i]],{i,Length[cs]}]]/.{\[CapitalDelta] r[1,a_]:>  1+ a r[1,a],\[CapitalDelta]^2 r[1,a_]:>  \[CapitalDelta]+a+ a^2 r[1,a]}/.{\[CapitalDelta] r[2,a_]:>  a r[2,a]+ r[1,a],\[CapitalDelta]^2 r[2,a_]:>  1+a^2 r[2,a]+2a r[1,a],\[CapitalDelta]^3 r[2,a_]:>2a+\[CapitalDelta]+a^3 r[2,a]+3a^2 r[1,a]};
,{n,1,nmax},{m,0,2(nmax-n)+mmax}];
];
Clear[z];
Return[];
);
];




Options[Recurse2]={coeffTab->None};
Recurse2[eps_,ll_,Delta12_,Delta34_,prec_,mmax_,nmax_,kmax_,OptionsPattern[]]:=Module[{\[Epsilon],l,fudge,expr,zdermax,t1,t2,t3,t4,coeffs,cn,cnRhoDer,CRec,cs,\[Rho]num,coeffTb,z,casCoeffs},(

ClearAll[cn,cnRhoDer,CRec];

\[Epsilon]=eps;l=ll;
fudge=10^-50;

coeffs[n_]:={(n^4+4 n^3 (-6+\[CapitalDelta]+\[Epsilon])-(6+l-\[CapitalDelta]) (-7+2 \[CapitalDelta]) (-5+2 \[Epsilon]) (-6+l+\[CapitalDelta]+2 \[Epsilon])-2 n (-6+\[CapitalDelta]+\[Epsilon]) (-71+l^2-(-22+\[CapitalDelta]) \[CapitalDelta]+2 (13+l-3 \[CapitalDelta]) \[Epsilon])+n^2 (-l^2+5 (43+(-14+\[CapitalDelta]) \[CapitalDelta])-2 (37+l-7 \[CapitalDelta]) \[Epsilon]+4 \[Epsilon]^2))/(n (-1+l+n+\[CapitalDelta]) (-1-l+n+\[CapitalDelta]-2 \[Epsilon]) (-2+n+2 \[CapitalDelta]-2 \[Epsilon])),(n^4+4 n^3 (-5-2 Delta12+2 Delta34+\[CapitalDelta]+\[Epsilon])+n^2 (149-132 Delta34-l^2-58 \[CapitalDelta]+5 \[CapitalDelta]^2+4 Delta12 (33+Delta34-6 \[CapitalDelta]-6 \[Epsilon])-62 \[Epsilon]-2 l \[Epsilon]+14 \[CapitalDelta] \[Epsilon]+4 \[Epsilon]^2+24 Delta34 (\[CapitalDelta]+\[Epsilon]))-4 ((5+l-\[CapitalDelta]) (-3+\[CapitalDelta]) (-2+\[Epsilon]) (-5+l+\[CapitalDelta]+2 \[Epsilon])+1/2 Delta34 (-11+2 \[CapitalDelta]+2 \[Epsilon]) (-60+l^2-(-20+\[CapitalDelta]) \[CapitalDelta]+2 (12+l-3 \[CapitalDelta]) \[Epsilon])-1/2 Delta12 (-11+2 \[CapitalDelta]+2 \[Epsilon]) (-60+l^2-(-20+\[CapitalDelta]) \[CapitalDelta]+2 (12+l-3 \[CapitalDelta]) \[Epsilon]+Delta34 (-6+\[CapitalDelta]+2 \[Epsilon])))-2 n ((-5+\[CapitalDelta]+\[Epsilon]) (-49+l^2-(-18+\[CapitalDelta]) \[CapitalDelta]+2 (11+l-3 \[CapitalDelta]) \[Epsilon])-2 Delta34 (181-l^2+\[CapitalDelta] (-64+5 \[CapitalDelta])-2 (34+l-7 \[CapitalDelta]) \[Epsilon]+4 \[Epsilon]^2)-2 Delta12 (-181+l^2+(64-5 \[CapitalDelta]) \[CapitalDelta]+2 (34+l-7 \[CapitalDelta]) \[Epsilon]-4 \[Epsilon]^2+1/2 Delta34 (-23+4 \[CapitalDelta]+6 \[Epsilon]))))/(n (-1+l+n+\[CapitalDelta]) (-1-l+n+\[CapitalDelta]-2 \[Epsilon]) (-2+n+2 \[CapitalDelta]-2 \[Epsilon])),(-1620+65 l^2-3 n^4+2 \[CapitalDelta]^3 (11-4 \[Epsilon])+4 n^3 (14-2 Delta12+2 Delta34-3 \[CapitalDelta]-\[Epsilon])+580 \[Epsilon]+130 l \[Epsilon]-18 l^2 \[Epsilon]+40 \[Epsilon]^2-36 l \[Epsilon]^2+\[CapitalDelta]^2 (-277+102 \[Epsilon])-4 Delta34^2 (l^2-(-8+\[CapitalDelta]) (-15+4 \[CapitalDelta])+2 (20+l-5 \[CapitalDelta]) \[Epsilon])+2 Delta34 (-9+2 \[CapitalDelta]+2 \[Epsilon]) (40-l^2+(-16+\[CapitalDelta]) \[CapitalDelta]-2 (10+l-3 \[CapitalDelta]) \[Epsilon])+2 \[CapitalDelta] (584-\[Epsilon] (219+2 \[Epsilon])+l^2 (-11+4 \[Epsilon])+2 l \[Epsilon] (-11+4 \[Epsilon]))-4 Delta12^2 (-120+l^2+47 \[CapitalDelta]-4 \[CapitalDelta]^2+2 (20+l-5 \[CapitalDelta]) \[Epsilon]+Delta34 (-21+4 \[CapitalDelta]+6 \[Epsilon]))+n^2 (-401+20 Delta12^2+20 Delta34^2+3 l^2-15 \[CapitalDelta]^2-18 \[CapitalDelta] (-9+\[Epsilon])+62 \[Epsilon]+6 l \[Epsilon]+4 \[Epsilon]^2+12 Delta34 (-9+2 \[CapitalDelta]+2 \[Epsilon])-4 Delta12 (-27+11 Delta34+6 \[CapitalDelta]+6 \[Epsilon]))-2 Delta12 (-2 Delta34^2 (-21+4 \[CapitalDelta]+6 \[Epsilon])+(-9+2 \[CapitalDelta]+2 \[Epsilon]) (40-l^2+(-16+\[CapitalDelta]) \[CapitalDelta]-2 (10+l-3 \[CapitalDelta]) \[Epsilon])+Delta34 (555-4 l^2-213 \[CapitalDelta]+18 \[CapitalDelta]^2+(-220-8 l+50 \[CapitalDelta]) \[Epsilon]+12 \[Epsilon]^2))-2 n (-242 Delta34+14 l^2+373 \[CapitalDelta]+3 \[CapitalDelta]^3+2 Delta12^2 (49+4 Delta34-10 \[CapitalDelta]-8 \[Epsilon])+163 (-4+\[Epsilon])+28 l \[Epsilon]-l^2 \[Epsilon]+14 \[Epsilon]^2-2 l \[Epsilon]^2-2 Delta34^2 (-49+10 \[CapitalDelta]+8 \[Epsilon])+\[CapitalDelta]^2 (-64+11 \[Epsilon])+2 Delta34 (l^2+(52-5 \[CapitalDelta]) \[CapitalDelta]+2 (28+l-7 \[CapitalDelta]) \[Epsilon]-4 \[Epsilon]^2)-2 Delta12 (-121+4 Delta34^2+l^2+(52-5 \[CapitalDelta]) \[CapitalDelta]+1/2 Delta34 (221-44 \[CapitalDelta]-42 \[Epsilon])+2 (28+l-7 \[CapitalDelta]) \[Epsilon]-4 \[Epsilon]^2)-\[CapitalDelta] (3 l^2+6 l \[Epsilon]+2 \[Epsilon] (44+\[Epsilon]))))/(n (-1+l+n+\[CapitalDelta]) (-1-l+n+\[CapitalDelta]-2 \[Epsilon]) (-2+n+2 \[CapitalDelta]-2 \[Epsilon])),(-3 n^4-4 n^3 (-11-4 Delta12+4 Delta34+3 \[CapitalDelta]+\[Epsilon])+n^2 (-251+20 Delta12^2+20 Delta34^2+3 l^2-48 Delta12 (4+Delta34-\[CapitalDelta])-48 Delta34 (-4+\[CapitalDelta])-15 \[CapitalDelta]^2-18 \[CapitalDelta] (-7+\[Epsilon])+50 \[Epsilon]+6 l \[Epsilon]+4 \[Epsilon]^2)+2 n (329-8 Delta12^3+8 Delta34^3-11 l^2-3 \[CapitalDelta]^3+\[CapitalDelta]^2 (49-11 \[Epsilon])-107 \[Epsilon]-22 l \[Epsilon]+l^2 \[Epsilon]-10 \[Epsilon]^2+2 l \[Epsilon]^2+2 Delta34^2 (-39+10 \[CapitalDelta]+8 \[Epsilon])+2 Delta12^2 (-39+16 Delta34+10 \[CapitalDelta]+8 \[Epsilon])+4 Delta34 (-99+l^2+(46-5 \[CapitalDelta]) \[CapitalDelta]+2 (l-\[CapitalDelta]) \[Epsilon]+4 \[Epsilon]^2)+\[CapitalDelta] (-229+3 l^2+6 l \[Epsilon]+2 \[Epsilon] (35+\[Epsilon]))-4 Delta12 (-99+8 Delta34^2+l^2+(46-5 \[CapitalDelta]) \[CapitalDelta]+2 (l-\[CapitalDelta]) \[Epsilon]+4 \[Epsilon]^2+1/2 Delta34 (-91+24 \[CapitalDelta]+14 \[Epsilon])))-4 (-4 Delta12^3 (4+Delta34-\[CapitalDelta])-4 Delta34^3 (-4+\[CapitalDelta])+\[CapitalDelta]^2 (41-20 \[Epsilon])+2 \[CapitalDelta]^3 (-2+\[Epsilon])+Delta34^2 (-76+l^2+(37-4 \[CapitalDelta]) \[CapitalDelta]+2 (16+l-5 \[CapitalDelta]) \[Epsilon])+\[CapitalDelta] (-143-2 l^2 (-2+\[Epsilon])+70 \[Epsilon]-4 l (-2+\[Epsilon]) \[Epsilon])+2 Delta34 (-4+\[CapitalDelta]) (35-l^2+\[CapitalDelta]^2+2 \[CapitalDelta] (-7+\[Epsilon])-2 l \[Epsilon]-4 \[Epsilon]^2)+2 (83-5 l^2+(-39+2 (-5+l) l) \[Epsilon]+(-2+4 l) \[Epsilon]^2)+Delta12^2 (-76+8 Delta34^2+l^2+(37-4 \[CapitalDelta]) \[CapitalDelta]+2 (16+l-5 \[CapitalDelta]) \[Epsilon]-Delta34 (-71+16 \[CapitalDelta]+10 \[Epsilon]))-2 Delta12 (2 Delta34^3-1/2 Delta34^2 (-71+16 \[CapitalDelta]+10 \[Epsilon])+(-4+\[CapitalDelta]) (35-l^2+\[CapitalDelta]^2+2 \[CapitalDelta] (-7+\[Epsilon])-2 l \[Epsilon]-4 \[Epsilon]^2)+1/2 Delta34 (-172+2 l^2+(87-10 \[CapitalDelta]) \[CapitalDelta]+2 (22+2 l-9 \[CapitalDelta]) \[Epsilon]+8 \[Epsilon]^2))))/(n (-1+l+n+\[CapitalDelta]) (-1-l+n+\[CapitalDelta]-2 \[Epsilon]) (-2+n+2 \[CapitalDelta]-2 \[Epsilon])),(468-33 l^2+3 n^4+16 Delta34^3 (-3+\[CapitalDelta])-16 Delta12^3 (-3+Delta34+\[CapitalDelta])+108 \[Epsilon]-66 l \[Epsilon]+2 l^2 \[Epsilon]-72 \[Epsilon]^2+4 l \[Epsilon]^2-4 n^3 (10-4 Delta12+4 Delta34-3 \[CapitalDelta]+\[Epsilon])+2 \[CapitalDelta]^3 (-7+2 \[Epsilon])+4 Delta34^2 (-48+l^2+(29-4 \[CapitalDelta]) \[CapitalDelta]+2 (l+3 (-4+\[CapitalDelta])) \[Epsilon])+4 Delta12^2 (-48+8 Delta34^2+l^2+(29-4 \[CapitalDelta]) \[CapitalDelta]+1/2 Delta34 (-82+32 \[CapitalDelta]-20 \[Epsilon])+2 (-12+l+3 \[CapitalDelta]) \[Epsilon])+8 Delta34 (-3+\[CapitalDelta]) (l^2-(-7+\[CapitalDelta]) (-3+\[CapitalDelta])+2 (l-\[CapitalDelta]) \[Epsilon]+4 \[Epsilon]^2)-2 \[CapitalDelta] (216-7 l^2+(11+2 (-7+l) l) \[Epsilon]+(-22+4 l) \[Epsilon]^2)+n^2 (209-20 Delta12^2-20 Delta34^2-3 l^2-48 Delta34 (-3+\[CapitalDelta])+15 \[CapitalDelta]^2+48 Delta12 (-3+Delta34+\[CapitalDelta])+34 \[Epsilon]-6 l \[Epsilon]-4 \[Epsilon]^2-6 \[CapitalDelta] (19+\[Epsilon]))+\[CapitalDelta]^2 (133-2 \[Epsilon] (7+4 \[Epsilon]))-8 Delta12 (2 Delta34^3+1/4 Delta34^2 (-82+32 \[CapitalDelta]-20 \[Epsilon])+(-3+\[CapitalDelta]) (l^2-(-7+\[CapitalDelta]) (-3+\[CapitalDelta])+2 (l-\[CapitalDelta]) \[Epsilon]+4 \[Epsilon]^2)+1/2 Delta34 (-123+2 l^2+(73-10 \[CapitalDelta]) \[CapitalDelta]+2 (-27+2 l+5 \[CapitalDelta]) \[Epsilon]+8 \[Epsilon]^2))+2 n (-252-8 Delta12^3+8 Delta34^3+10 l^2+3 \[CapitalDelta]^3+\[CapitalDelta]^2 (-44+\[Epsilon])-51 \[Epsilon]+20 l \[Epsilon]+l^2 \[Epsilon]+18 \[Epsilon]^2+2 l \[Epsilon]^2+2 Delta34^2 (31-10 \[CapitalDelta]+8 \[Epsilon])+2 Delta12^2 (31+16 Delta34-10 \[CapitalDelta]+8 \[Epsilon])-3 \[CapitalDelta] (-63+l^2+2 l \[Epsilon]+2 (-2+\[Epsilon]) \[Epsilon])+4 Delta34 (-57+l^2+34 \[CapitalDelta]-5 \[CapitalDelta]^2+2 l \[Epsilon]-2 \[CapitalDelta] \[Epsilon]+4 \[Epsilon]^2)-4 Delta12 (-57+(77 Delta34)/2+8 Delta34^2+l^2-5 \[CapitalDelta]^2+7 Delta34 \[Epsilon]+2 l \[Epsilon]+4 \[Epsilon]^2-2 \[CapitalDelta] (-17+6 Delta34+\[Epsilon]))))/(n (-1+l+n+\[CapitalDelta]) (-1-l+n+\[CapitalDelta]-2 \[Epsilon]) (-2+n+2 \[CapitalDelta]-2 \[Epsilon])),(3 n^4+4 n^3 (-7-2 Delta12+2 Delta34+3 \[CapitalDelta]-\[Epsilon])+n^2 (107-20 Delta12^2-20 Delta34^2-3 l^2+15 \[CapitalDelta]^2+22 \[Epsilon]-6 l \[Epsilon]-4 \[Epsilon]^2-6 \[CapitalDelta] (13+\[Epsilon])-12 Delta34 (5-2 \[CapitalDelta]+2 \[Epsilon])+4 Delta12 (15+11 Delta34-6 \[CapitalDelta]+6 \[Epsilon]))-2 n (97-7 l^2-3 \[CapitalDelta]^3+2 Delta12^2 (-21+4 Delta34+10 \[CapitalDelta]-8 \[Epsilon])-\[CapitalDelta]^2 (-29+\[Epsilon])+23 \[Epsilon]-14 l \[Epsilon]-l^2 \[Epsilon]-14 \[Epsilon]^2-2 l \[Epsilon]^2-2 Delta34^2 (21-10 \[CapitalDelta]+8 \[Epsilon])+3 \[CapitalDelta] (-31+l^2+2 l \[Epsilon]+2 (-1+\[Epsilon]) \[Epsilon])-2 Delta12 (-37+4 Delta34^2+l^2+(28-5 \[CapitalDelta]) \[CapitalDelta]+1/2 Delta34 (-87+44 \[CapitalDelta]-42 \[Epsilon])+2 (-14+l+5 \[CapitalDelta]) \[Epsilon]-4 \[Epsilon]^2)-2 Delta34 (37-l^2+\[CapitalDelta] (-28+5 \[CapitalDelta])-2 (-14+l+5 \[CapitalDelta]) \[Epsilon]+4 \[Epsilon]^2))+4 (34-4 l^2+\[CapitalDelta]^3 (-2+\[Epsilon])+9 \[Epsilon]-8 l \[Epsilon]+l^2 \[Epsilon]-10 \[Epsilon]^2+2 l \[Epsilon]^2+1/2 Delta34 (5-2 \[CapitalDelta]+2 \[Epsilon]) (l^2-(-6+\[CapitalDelta]) (-2+\[CapitalDelta])+2 (-4+l+\[CapitalDelta]) \[Epsilon])+Delta34^2 (-22+l^2+19 \[CapitalDelta]-4 \[CapitalDelta]^2+2 (-8+l+3 \[CapitalDelta]) \[Epsilon])-\[CapitalDelta] (39+l^2 (-2+\[Epsilon])+\[Epsilon]+2 l (-2+\[Epsilon]) \[Epsilon]-8 \[Epsilon]^2)+\[CapitalDelta]^2 (15-\[Epsilon] (3+2 \[Epsilon]))+Delta12^2 (-22+l^2+19 \[CapitalDelta]-4 \[CapitalDelta]^2+2 (-8+l+3 \[CapitalDelta]) \[Epsilon]+Delta34 (7-4 \[CapitalDelta]+6 \[Epsilon]))-1/2 Delta12 (2 Delta34^2 (7-4 \[CapitalDelta]+6 \[Epsilon])+(5-2 \[CapitalDelta]+2 \[Epsilon]) (l^2-(-6+\[CapitalDelta]) (-2+\[CapitalDelta])+2 (-4+l+\[CapitalDelta]) \[Epsilon])-Delta34 (86-4 l^2-79 \[CapitalDelta]+18 \[CapitalDelta]^2+(74-8 l-34 \[CapitalDelta]) \[Epsilon]+12 \[Epsilon]^2))))/(n (-1+l+n+\[CapitalDelta]) (-1-l+n+\[CapitalDelta]-2 \[Epsilon]) (-2+n+2 \[CapitalDelta]-2 \[Epsilon])),(-n^4+4 n^3 (2-2 Delta12+2 Delta34-\[CapitalDelta]+\[Epsilon])+n^2 (-23+l^2+22 \[CapitalDelta]-5 \[CapitalDelta]^2-4 Delta12 (-9+Delta34+6 \[CapitalDelta]-6 \[Epsilon])-22 \[Epsilon]+2 l \[Epsilon]+10 \[CapitalDelta] \[Epsilon]-4 \[Epsilon]^2-12 Delta34 (3-2 \[CapitalDelta]+2 \[Epsilon]))+(3-2 \[CapitalDelta]+2 \[Epsilon]) ((1+2 Delta34) (-2+l+\[CapitalDelta]) (2+l-\[CapitalDelta]+2 \[Epsilon])-2 Delta12 ((-2+l+\[CapitalDelta]) (2+l-\[CapitalDelta]+2 \[Epsilon])+1/2 Delta34 (2-2 \[CapitalDelta]+4 \[Epsilon])))-2 n ((2-\[CapitalDelta]+\[Epsilon]) (-7+l^2-(-6+\[CapitalDelta]) \[CapitalDelta]+2 (-3+l+\[CapitalDelta]) \[Epsilon])+2 Delta34 (-13+l^2+(16-5 \[CapitalDelta]) \[CapitalDelta]+2 (-8+l+5 \[CapitalDelta]) \[Epsilon]-4 \[Epsilon]^2)-2 Delta12 (-13+l^2+(16-5 \[CapitalDelta]) \[CapitalDelta]+2 (-8+l+5 \[CapitalDelta]) \[Epsilon]-4 \[Epsilon]^2+1/2 Delta34 (5-4 \[CapitalDelta]+6 \[Epsilon]))))/(n (-1+l+n+\[CapitalDelta]) (-1-l+n+\[CapitalDelta]-2 \[Epsilon]) (-2+n+2 \[CapitalDelta]-2 \[Epsilon]))};

zdermax=2nmax+mmax;


\[Rho]num=SetPrecision[3-2Sqrt[2],prec];

z=1/2; (*If we want to generate blocks at some other point this needs to be changed*)
coeffTb=OptionValue[coeffTab];

casCoeffs[m_,n_,\[Epsilon]_,l_]={-(((-2+m) (-1+m) m)/(8 (-1+z) z^2)),-(((-1+m) m (-1+3 z))/(4 (-1+z) z^2)),-((m (22-13 m+m^2-34 n+12 m n+12 n^2+Delta34 (-5+m+4 n)-Delta12 (-5+Delta34+m+4 n)+2 (1+m-2 n) \[Epsilon]))/(16 (-1+z) z^2 (-1+2 n+2 \[Epsilon]))),(m (2-3 z))/(2 (-1+z) z),-(1/(8 (-1+z) z^2 (-1+2 n+2 \[Epsilon])))(-6+l^2-(-10+Delta12 (-4+Delta34)+4 Delta34) z+m^2 (-1+3 z)+2 n (5-2 n+(-11-2 Delta12+2 Delta34+6 n) z)+m (9-8 n+(-23-2 Delta12+2 Delta34+24 n) z)+(-2+\[CapitalDelta]) \[CapitalDelta]+2 (-2+l-2 n (-1+z)+2 (1+m) z-\[CapitalDelta]) \[Epsilon]),((-1+n) (6+Delta12-Delta34-3 m-4 n+2 \[Epsilon]))/(16 (-1+z) z^2 (-1+2 n+2 \[Epsilon])),(-8+8 n+m (2-3 z)+(10+Delta12-Delta34) z-2 z (6 n+\[Epsilon]))/(4 (-1+z) z (-1+2 n+2 \[Epsilon])),-(((-1+n) (-1+3 z))/(8 (-1+z) z^2 (-1+2 n+2 \[Epsilon]))),1/(2-4 n-4 \[Epsilon])}//Simplify;

casInds[m_,n_]={{-3+m,n},{-2+m,n},{-1+m,-1+n},{-1+m,n}, {m,-1+n},{1+m,-2+n},{1+m,-1+n}, {2+m,-2+n},{2+m,-1+n}};

(*************************)





cn[n_]:=cn[n]=Chop[Apart[Expand[Sum[coeffs[n][[i]]cn[n-8+i],{i,7}]],\[CapitalDelta]]];
cn[0]=1;
cn[x_/;x<0]=0;


(* Power series are rho^\Delta Sum( cnrhoder[n,i] rho^(n-i));  *)
(* fudge factor is required to prevent certain cancellations from occurring; otherwise different components will have different poles *)



cnRhoDer[n_,i_]:=cnRhoDer[n,i]=Expand[(\[CapitalDelta]+n-(i-1)+fudge)cnRhoDer[n,i-1]]/.{\[CapitalDelta] r[1,a_]:>  1+ a r[1,a],\[CapitalDelta]^2 r[1,a_]:>  \[CapitalDelta]+a+ a^2 r[1,a]}/.{\[CapitalDelta] r[2,a_]:>  a r[2,a]+ r[1,a],\[CapitalDelta]^2 r[2,a_]:>  1+a^2 r[2,a]+2a r[1,a],\[CapitalDelta]^3 r[2,a_]->2a+\[CapitalDelta]+a^3 r[2,a]+3a^2 r[1,a]};
cnRhoDer[n_,0]:=cnRhoDer[n,0]=cn[n]/.{1/\[CapitalDelta]:>r[1,0],1/\[CapitalDelta]^2:>r[2,0],1/(a_+b_ \[CapitalDelta]):> 1/b r[1,-a/b],1/(a_+\[CapitalDelta]):> r[1,-a]}/.{Power[a_+b_ \[CapitalDelta],m_]:> 1/b^-m r[-m,-a/b],Power[a_+\[CapitalDelta],m_]:> r[-m,-a]};

t1=AbsoluteTiming[
Do[cn[k],{k,0,kmax}]];


(*Need to convert to z=zb derivatives (or more precisely the a derivatives *)
t2=AbsoluteTiming[
AllDer[ll,0,0]=(Sum[Expand[cn[k]\[Rho]num^k],{k,0,kmax}]//Expand)/.{1/\[CapitalDelta]:>r[1,0],1/\[CapitalDelta]^2:>r[2,0],1/(a_+b_ \[CapitalDelta]):> 1/b r[1,-a/b],1/(a_+\[CapitalDelta]):> r[1,-a]}/.{Power[a_+b_ \[CapitalDelta],m_]:> 1/b^-m r[-m,-a/b],Power[a_+\[CapitalDelta],m_]:> r[-m,-a]}//Simplify;
];

t3=AbsoluteTiming[
Do[
AllDer[ll,m,0]=(1/2^m Sum[Expand[coeffTb[[m,i]]cnRhoDer[k,i]\[Rho]num^(k-i)]/.{\[CapitalDelta] r[1,a_]:>  1+ a r[1,a],\[CapitalDelta]^2 r[1,a_]:>  \[CapitalDelta]+a+ a^2 r[1,a]}/.{\[CapitalDelta] r[2,a_]:>  a r[2,a]+ r[1,a],\[CapitalDelta]^2 r[2,a_]:>  1+a^2 r[2,a]+2a r[1,a],\[CapitalDelta]^3 r[2,a_]->2a+\[CapitalDelta]+a^3 r[2,a]+3a^2 r[1,a]},{k,0,kmax},{i,1,Length[coeffTb[[m]]]}]);
,{m,1,zdermax}];
];

(* Use Casimir rec*)
CRec[m_/;m<0,n_]:=0;
CRec[m_,n_/;n<0]:=0;
Do[CRec[m,0]=AllDer[ll,m,0],{m,0,zdermax}];

t4=AbsoluteTiming[
Do[
cs=casCoeffs[m,n,eps,ll]//Simplify;

AllDer[ll,m,n]=CRec[m,n]=Expand[Sum[cs[[i]]CRec@@casInds[m,n][[i]],{i,Length[cs]}]]/.{\[CapitalDelta] r[1,a_]:>  1+ a r[1,a],\[CapitalDelta]^2 r[1,a_]:>  \[CapitalDelta]+a+ a^2 r[1,a]}/.{\[CapitalDelta] r[2,a_]:>  a r[2,a]+ r[1,a],\[CapitalDelta]^2 r[2,a_]:>  1+a^2 r[2,a]+2a r[1,a],\[CapitalDelta]^3 r[2,a_]:>2a+\[CapitalDelta]+a^3 r[2,a]+3a^2 r[1,a]};
,{n,1,nmax},{m,0,2(nmax-n)+mmax}];
];
Clear[z];
Return[];
);
];




GetDataN0[L_,\[Epsilon]_,nmax_,mmax_]:=Module[{tmp,Poles1,PoleCoeffs1,Poles2,PoleCoeffs2,cases1,cases2,Poly},(
Do[
tmp=2^(m+2n) AllDer[L,m,n]/.r[1,x_]->r1[x]/.r[2,x_]->r2[x]//Expand; (* The powers of 2 are to match Sheer and Slava's conventions. *)
(*Poly[m,n]=CoefficientList[tmp/.r1[x_]->0/.r2[x_]->0/.\[CapitalDelta]->\[CapitalDelta]+2\[Epsilon]+L,\[CapitalDelta]]/.a_/;MantissaExponent[a][[1]]==0.->0;*)

Poly[m,n]=CoefficientList[tmp/.r1[x_]->0/.r2[x_]->0,\[CapitalDelta]]/.a_/;MantissaExponent[a][[1]]==0.->0;


cases1=Cases[tmp,a_ r1[x_]];

Poles1[m,n]=cases1/.a_ r1[x_]->x;
PoleCoeffs1[m,n]=cases1/.a_ r1[x_]->a;

cases2=Cases[tmp,a_ r2[x_]];

Poles2[m,n]=cases2/.a_ r2[x_]->x;
PoleCoeffs2[m,n]=cases2/.a_ r2[x_]->a;



,{n,0,nmax},{m,0,2(nmax-n)+mmax}];

Return[{Poly,Poles1,PoleCoeffs1,Poles2,PoleCoeffs2}];

);
];


Options[ComputeTable]={kMax->60,Prec->64,AllSpins->False,\[CapitalDelta]12->0,\[CapitalDelta]34->0};
ComputeTable[Dim_,nmax_,mmax_,Lmax_,filename_,OptionsPattern[]]:=Module[{allspins,spinstep,eps,Lct,Res,t1,t2,t3,t4,t5,polylist,polelist1,polelist2,clist1,clist2,plen,LL,maxPolylength,L,\[CapitalDelta]0,Poly,Poles1,PoleCoeffs1,Poles2,PoleCoeffs2,file,prec,kmax,Delta12,Delta34,zdermax,z,coeffTb},

prec=OptionValue[Prec];
kmax=OptionValue[kMax];
allspins=OptionValue[AllSpins];
spinstep=If[allspins,1,2];
eps=Rationalize[(Dim-2)/2];
Delta12=OptionValue[\[CapitalDelta]12];
       Delta34=OptionValue[\[CapitalDelta]34];

zdermax=2nmax+mmax;

z=1/2; (*If we want to generate blocks at some other point this needs to be changed*)
coeffTb=SetNDer[zdermax,prec,z];


If[(Delta12!=0||Delta34!=0),
Monitor[For[Lct=0,Lct<=Lmax,Lct+=spinstep,
Recurse2[eps,Lct,Delta12,Delta34,prec,mmax,nmax,kmax,coeffTab->coeffTb];]
,Lct];
,
Monitor[For[Lct=0,Lct<=Lmax,Lct+=spinstep,
Recurse[eps,Lct,prec,mmax,nmax,kmax,coeffTab->coeffTb];]
,Lct];
];

Res={};For[LL=0,LL<=Lmax,LL+=spinstep,
{Poly,Poles1,PoleCoeffs1,Poles2,PoleCoeffs2}=GetDataN0[LL,eps,nmax,mmax];
t1=Flatten[Table[Poly[m,n],{n,0,nmax},{m,0,2(nmax-n)+mmax}],1];
t2=DeleteDuplicates[Flatten[Table[Poles1[0,0],{n,0,nmax},{m,0,2(nmax-n)+mmax}]]];
t3=Flatten[Table[PoleCoeffs1[m,n],{n,0,nmax},{m,0,2(nmax-n)+mmax}],1];
t4=DeleteDuplicates[Flatten[Table[Poles2[0,0],{n,0,nmax},{m,0,2(nmax-n)+mmax}]]];
t5=Flatten[Table[PoleCoeffs2[m,n],{n,0,nmax},{m,0,2(nmax-n)+mmax}],1];


AppendTo[Res,{t1,t2,t3,t4,t5}];
];

plen=Length[Flatten[Table[0, {n,0,nmax},{m,0,mmax+2(nmax-n)}]]];
file=OpenWrite[filename];
Write[file,N[eps]];
Write[file, nmax];
Write[file, mmax];
Write[file,Lmax];
Write[file,Boole[allspins]]; (* whether there are odd spins; 1 means yes, 0 means no *)
Write[file,Round[prec Log[2,10]]-1]; (*nr of binary digits of precision *)
Write[file,plen];

For[Lct=1,Lct<=Length[Res],Lct+=1,

polylist=Res[[Lct,1]];
polelist1=Res[[Lct,2]];
clist1=Res[[Lct,3]];
polelist2=Res[[Lct,4]];
clist2=Res[[Lct,5]];


plen=Length[Flatten[Table[0, {n,0,nmax},{m,0,mmax+2(nmax-n)}]]];
maxPolylength=Max[Table[Length[polylist[[i]]],{i,plen}]];
L=spinstep*(Lct-1);
\[CapitalDelta]0=2eps+L;

WriteString[file,ToString[CForm[N[\[CapitalDelta]0,prec]]]<>"\n"];
Write[file,maxPolylength];

Do[Map[WriteString[file,(ToString[CForm[#]]<>"\n")]&,PadRight[polylist[[i]],maxPolylength,0]];,{i,plen}];

Write[file,Length[polelist1]];
Map[WriteString[file,(ToString[CForm[#]]<>"\n")]&,N[polelist1,prec]];
Do[Map[WriteString[file,(ToString[CForm[#]]<>"\n")]&,clist1[[i]]];,{i,plen}];

Write[file,Length[polelist2]];
Map[WriteString[file,(ToString[CForm[#]]<>"\n")]&,N[polelist2,prec]];
Do[Map[WriteString[file,(ToString[CForm[#]]<>"\n")]&,clist2[[i]]];,{i,plen}];

];


Close[file];

];

