
#include <sys/types.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "resource.h"

#include "f2c.h"
#include <string.h>
//#include <stdlib.h>
#pragma comment(lib, "vcf2c.lib")

/* Table of constant values */

static integer c__1 = 1;
static integer c__3 = 3;
static integer c__9 = 9;
static integer c__4 = 4;

//主函数，读取输入文件和数据文件输出总剂量
__declspec(dllexport) double __stdcall dose(void)
{
    /* Initialized data */

    static char det[8 * 11] = "Aluminum" "Graphite" "Silicon " "Air     " "Bon"
        "e    " "CaF2    " "GaAs    " "LiF     " "SiO2    " "Tissue  "
        "H2O     ";
    static real zmin = 1e-6f;
    static real radcon = 1.6021892e-8f;
    static integer nbege = 1;
    static real enmu = .03f;
    static char version[4] = "2.10";

    /* Format strings */
    static char fmt_10[] = "(\002 Enter input filename: \002)";
    static char fmt_20[] = "(a)";
    static char fmt_30[] = "(\002 OUTPUT FROM SHIELDOSE-2, VERSION \002,a)";
    static char fmt_32[] = "(\002 Input filename: \002,a)";
    static char fmt_34[] = "(\002 Print-out filename: \002,a)";
    static char fmt_36[] = "(\002 Array output filename: \002,a)";
    static char fmt_370[] = "(/\002  IDET  INUC  IMAX  IUNT\002)";
    static char fmt_380[] = "(12i6)";
    static char fmt_450[] = "(/\002 SHIELD DEPTH (mils)\002)";
    static char fmt_455[] = "(1p6e12.5)";
    static char fmt_480[] = "(/\002 SHIELD DEPTH (g/cm2)\002)";
    static char fmt_510[] = "(/\002 SHIELD DEPTH (mm)\002)";
    static char fmt_1435[] = "(1p10e10.3)";
    static char fmt_550[] = "(/\002     EMINS     EMAXS     EMINP     EMAXP "
        "NPTSP     EMINE     EMAXE NPTSE\002)";
    static char fmt_560[] = "(4f10.3,i6,2f10.3,i6)";
    static char fmt_580[] = "(4f10.3,i6,2f10.3,i6,\002  ADJUSTED VALUES\002)";
    static char fmt_840[] = "(/)";
    static char fmt_850[] = "(\002 \002)";
    static char fmt_860[] = "(4x,a)";
    static char fmt_870[] = "(/\002 JSMAX JPMAX JEMAX       EUNIT      DUR"
        "ATN\002)";
    static char fmt_880[] = "(3i6,1p2e12.5)";
    static char fmt_885[] = "(//\002 E(MeV)\002)";
    static char fmt_905[] = "(1p10e12.4)";
    static char fmt_890[] = "(/\002 SOLAR PROTON SPECTRUM (/energy/cm2)\002)";
    static char fmt_891[] = "(/\002 SPECTRUM INTEGRATED FROM\002,1pe11.4,"
        "\002 TO\002,1pe11.4,\002 MeV, USING\002,i5,\002 POINTS\002)";
    static char fmt_896[] = "(/\002 ASSUMED FRACTION OF BEAM ENERGY INTO NEU"
        "TRON ENERGY =\002,1pe12.5)";
    static char fmt_910[] = "(/\002 TRAPPED PROTON SPECTRUM (/energy/cm2/tim"
        "e)\002)";
    static char fmt_930[] = "(/\002 ELECTRON SPECTRUM (/energy/cm2/time)\002)"
        ;
    static char fmt_1180[] = "(//\002 DOSE AT TRANSMISSION SURFACE OF FINITE"
        " ALUMINUM SLAB SHIELDS\002)";
    static char fmt_1200[] = "(//\002 DOSE IN SEMI-INFINITE ALUMINUM MEDIU"
        "M\002)";
    static char fmt_1230[] = "(/\002 rads \002,a)";
    static char fmt_1240[] = "(/\002 Proton results without nuclear attenuat"
        "ion\002)";
    static char fmt_1250[] = "(/\002 Proton results with approximate treatme"
        "nt of nuclear attenuation\002)";
    static char fmt_1260[] = "(\002    neglecting transport of energy by neu"
        "trons\002)";
    static char fmt_1270[] = "(\002    and crude exponential transport of en"
        "ergy by neutrons\002)";
    static char fmt_1310[] = "(/\002    Z(mils)      Z(mm)   Z(g/cm2)   ELEC"
        "TRON    BREMS      EL+BR     TRP PROT   SOL PROT  EL+BR+TRP    T"
        "OTAL\002)";
    static char fmt_1320[] = "(1p10e11.3)";
    static char fmt_1350[] = "(//\002 1/2 DOSE AT CENTER OF ALUMINUM SPHERE"
        "S\002)";

    /* System generated locals */
    integer i__1, i__2, i__3;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer s_wsfe(cilist*), e_wsfe(void), s_rsfe(cilist*), do_fio(integer*
        , char*, ftnlen), e_rsfe(void), f_open(olist*), s_rsle(cilist*)
        , do_lio(integer*, integer*, char*, ftnlen), e_rsle(void),
        s_wsle(cilist*), e_wsle(void), f_clos(cllist*);
    double log(doublereal), exp(doublereal);
    //    /* Subroutine */ int s_stop(char *, ftnlen);

        /* Local variables */
    static real g[1001];
    static integer i__, j, k, l, m, n;
    static real s[301], z__[71], gb[142142]	/* was [1001][71][2] */, ee[
        81], ge[142142]	/* was [1001][71][2] */, ar[14];
        static integer ne;
        static real ep[133], bs[34], re[81], gp[71071]	/* was [1001][71] */,
            te[14], zb[47], ye[81];
        static integer np;
        static real rp[133], tp[49], rs[14], zl[71], zm[71], zs[37], dee, arb[14],
            arc[14], ard[14], reb[81], rec[81], red[81];
        static char tag[72];
        static real tee[1001], din[51], yeb[81], yec[81], yed[81], seg[1001], rpb[
            133], rpc[133], rpd[133], dum[51], rsb[14], rsc[14], rsd[14], tel[
                1001], eps[301], sol[1001], spg[1001], tpl[1001], zre[51], zmm[71]
                    , tpp[1001], dep, ans, eav, zrp[51], dalb[518]	/* was [14][
                    37] */, dale[476]	/* was [14][34] */, dele, dinb[51], dinc[51],
                    dind[51];
                static integer ilec;
                static real dalp[2499]	/* was [49][51] */;
                extern /* Subroutine */ int lcof_(integer*, real*, real*, real*, real
                    *, real*);
                static real delp;
                static integer idet;
                static real dosb[284]	/* was [71][2][2] */, fepn[30], ares[1001],
                    dose[284]	/* was [71][2][2] */, drin[51], rine[1001], ylde[1001]
                    ;
                static integer inuc, imax;
                extern /* Subroutine */ int scof_(integer*, real*, real*, real*, real
                    *, real*);
                static real ansr, dosp[142]	/* was [71][2] */, tepn[30];
                static integer imix;
                static real rinp;
                static integer isol;
                static real rins[1001], dost;
                static integer itrp, iunt;
                static real zrin, dalbb[518]	/* was [14][37] */, dalbc[518]	/*
                    was [14][37] */, dalbd[518]	/* was [14][37] */, daleb[476]	/*
                    was [14][34] */, dalec[476]	/* was [14][34] */, daled[476]	/*
                    was [14][34] */, dalpb[2499]	/* was [49][51] */, dalpc[
                        2499]	/* was [49][51] */, dalpd[2499]	/* was [49][51] */,
                            fepnb[30], fepnc[30], fepnd[30], dratb[1316]	/* was [14][
                            47][2] */, emine, drate[1428]	/* was [14][51][2] */, emaxe;
                        static integer nlene;
                        static real doseb;
                        static integer lmaxb, jemax, lmaxe, mmaxe;
                        static real dratp[2499]	/* was [49][51] */;
                        static integer nmaxe;
                        static real emins, emaxs, eminp, dosol[142]	/* was [71][2] */;
                        static integer inatt;
                        static real enewt[1001];
                        static integer mmaxp, kmaxp, nmaxp, inewt, lmaxp, lmaxs, lmaxt;
                        static real emaxp;
                        static integer nptse;
                        static real eminu, emaxu;
                        extern /* Subroutine */ int bspol_(real*, integer*, real*, real*,
                            real*, real*, real*, real*);
                        static real benmu, zrinl;
                        static integer jsmax, nptsp, jpmax;
                        static real eunit;
                        static integer nlens, nfsts, nlsts;
                        extern /* Subroutine */ int integ_(real*, real*, integer*, real*);
                        static real eneut;
                        static integer nlenp, nfstp, nlstp, nfste, nlste;
                        static real dratbb[1316]	/* was [14][47][2] */, dratbc[1316]	/*
                            was [14][47][2] */, dratbd[1316]	/* was [14][47][2] */, drateb[
                                1428]	/* was [14][51][2] */, dratec[1428]	/* was [14][
                                51][2] */, drated[1428]	/* was [14][51][2] */, deltae, eminel;
                                static char filenm[40];
                                static real deltap, dratpb[2499]	/* was [49][51] */, dratpc[2499]
                                    /* was [49][51] */;
                                static char arrfil[40];
                                static real dratpd[2499]	/* was [49][51] */;
                                extern /* Subroutine */ int eindex_(real*, real*, integer*, real*,
                                    real*, integer*, integer*, integer*);
                                static integer nlensb, nlenpb;
                                static real deltas;
                                extern /* Subroutine */ int sphere_(real*, real*, integer*, real*);
                                static real eminul;
                                static char prtfil[40];
                                static integer nfstsb, nfstpb, nlstpb, nlstsb;
                                static real duratn;
                                extern /* Subroutine */ int spectr_(integer*, real*, real*, real*,
                                    real*, real*, integer*, integer*, integer*, integer*, real*
                                    , real*, real*);
                                static real dosebp;

                                /* Fortran I/O blocks */
                                static cilist io___7 = { 0, 6, 0, fmt_10, 0 };
                                static cilist io___8 = { 0, 5, 0, fmt_20, 0 };
                                static cilist io___10 = { 0, 9, 0, fmt_20, 0 };
                                static cilist io___12 = { 0, 9, 0, fmt_20, 0 };
                                static cilist io___14 = { 0, 10, 0, fmt_30, 0 };
                                static cilist io___15 = { 0, 12, 0, fmt_30, 0 };
                                static cilist io___16 = { 0, 10, 0, fmt_32, 0 };
                                static cilist io___17 = { 0, 12, 0, fmt_32, 0 };
                                static cilist io___18 = { 0, 10, 0, fmt_34, 0 };
                                static cilist io___19 = { 0, 12, 0, fmt_34, 0 };
                                static cilist io___20 = { 0, 10, 0, fmt_36, 0 };
                                static cilist io___21 = { 0, 12, 0, fmt_36, 0 };
                                static cilist io___22 = { 0, 10, 0, fmt_370, 0 };
                                static cilist io___23 = { 0, 9, 0, 0, 0 };
                                static cilist io___28 = { 0, 10, 0, fmt_380, 0 };
                                static cilist io___29 = { 0, 12, 0, fmt_380, 0 };
                                static cilist io___32 = { 0, 6, 0, 0, 0 };
                                static cilist io___33 = { 0, 6, 0, 0, 0 };
                                static cilist io___34 = { 0, 11, 0, fmt_20, 0 };
                                static cilist io___36 = { 0, 11, 0, 0, 0 };
                                static cilist io___42 = { 0, 11, 0, 0, 0 };
                                static cilist io___45 = { 0, 11, 0, 0, 0 };
                                static cilist io___47 = { 0, 11, 0, 0, 0 };
                                static cilist io___50 = { 0, 11, 0, 0, 0 };
                                static cilist io___52 = { 0, 11, 0, 0, 0 };
                                static cilist io___55 = { 0, 11, 0, 0, 0 };
                                static cilist io___59 = { 0, 11, 0, 0, 0 };
                                static cilist io___62 = { 0, 11, 0, 0, 0 };
                                static cilist io___73 = { 0, 6, 0, 0, 0 };
                                static cilist io___74 = { 0, 11, 0, fmt_20, 0 };
                                static cilist io___75 = { 0, 11, 0, 0, 0 };
                                static cilist io___83 = { 0, 11, 0, 0, 0 };
                                static cilist io___85 = { 0, 11, 0, 0, 0 };
                                static cilist io___87 = { 0, 11, 0, 0, 0 };
                                static cilist io___89 = { 0, 11, 0, 0, 0 };
                                static cilist io___91 = { 0, 11, 0, 0, 0 };
                                static cilist io___93 = { 0, 11, 0, 0, 0 };
                                static cilist io___95 = { 0, 11, 0, 0, 0 };
                                static cilist io___97 = { 0, 11, 0, 0, 0 };
                                static cilist io___99 = { 0, 11, 0, 0, 0 };
                                static cilist io___101 = { 0, 11, 0, 0, 0 };
                                static cilist io___103 = { 0, 11, 0, 0, 0 };
                                static cilist io___105 = { 0, 11, 0, 0, 0 };
                                static cilist io___107 = { 0, 11, 0, 0, 0 };
                                static cilist io___109 = { 0, 11, 0, 0, 0 };
                                static cilist io___138 = { 0, 10, 0, fmt_450, 0 };
                                static cilist io___139 = { 0, 9, 0, 0, 0 };
                                static cilist io___141 = { 0, 10, 0, fmt_455, 0 };
                                static cilist io___144 = { 0, 10, 0, fmt_480, 0 };
                                static cilist io___145 = { 0, 9, 0, 0, 0 };
                                static cilist io___146 = { 0, 10, 0, fmt_455, 0 };
                                static cilist io___147 = { 0, 10, 0, fmt_510, 0 };
                                static cilist io___148 = { 0, 9, 0, 0, 0 };
                                static cilist io___149 = { 0, 10, 0, fmt_455, 0 };
                                static cilist io___151 = { 0, 12, 0, fmt_1435, 0 };
                                static cilist io___152 = { 0, 10, 0, fmt_550, 0 };
                                static cilist io___153 = { 0, 9, 0, 0, 0 };
                                static cilist io___162 = { 0, 10, 0, fmt_560, 0 };
                                static cilist io___177 = { 0, 10, 0, fmt_580, 0 };
                                static cilist io___178 = { 0, 12, 0, fmt_560, 0 };
                                static cilist io___179 = { 0, 6, 0, 0, 0 };
                                static cilist io___180 = { 0, 6, 0, 0, 0 };
                                static cilist io___192 = { 0, 6, 0, 0, 0 };
                                static cilist io___207 = { 0, 6, 0, 0, 0 };
                                static cilist io___208 = { 0, 10, 0, fmt_840, 0 };
                                static cilist io___209 = { 0, 10, 0, fmt_850, 0 };
                                static cilist io___210 = { 0, 9, 1, fmt_20, 0 };
                                static cilist io___211 = { 0, 6, 0, fmt_860, 0 };
                                static cilist io___212 = { 0, 10, 0, fmt_20, 0 };
                                static cilist io___213 = { 0, 12, 0, fmt_20, 0 };
                                static cilist io___214 = { 0, 10, 0, fmt_870, 0 };
                                static cilist io___215 = { 0, 9, 0, 0, 0 };
                                static cilist io___221 = { 0, 10, 0, fmt_880, 0 };
                                static cilist io___222 = { 0, 12, 0, fmt_880, 0 };
                                static cilist io___227 = { 0, 10, 0, fmt_885, 0 };
                                static cilist io___228 = { 0, 9, 0, 0, 0 };
                                static cilist io___231 = { 0, 10, 0, fmt_905, 0 };
                                static cilist io___232 = { 0, 12, 0, fmt_905, 0 };
                                static cilist io___233 = { 0, 10, 0, fmt_890, 0 };
                                static cilist io___234 = { 0, 9, 0, 0, 0 };
                                static cilist io___236 = { 0, 10, 0, fmt_905, 0 };
                                static cilist io___237 = { 0, 12, 0, fmt_905, 0 };
                                static cilist io___242 = { 0, 10, 0, fmt_891, 0 };
                                static cilist io___246 = { 0, 10, 0, fmt_896, 0 };
                                static cilist io___248 = { 0, 10, 0, fmt_885, 0 };
                                static cilist io___249 = { 0, 9, 0, 0, 0 };
                                static cilist io___250 = { 0, 10, 0, fmt_905, 0 };
                                static cilist io___251 = { 0, 12, 0, fmt_905, 0 };
                                static cilist io___252 = { 0, 10, 0, fmt_910, 0 };
                                static cilist io___253 = { 0, 9, 0, 0, 0 };
                                static cilist io___254 = { 0, 10, 0, fmt_905, 0 };
                                static cilist io___255 = { 0, 12, 0, fmt_905, 0 };
                                static cilist io___260 = { 0, 10, 0, fmt_891, 0 };
                                static cilist io___261 = { 0, 10, 0, fmt_896, 0 };
                                static cilist io___263 = { 0, 10, 0, fmt_885, 0 };
                                static cilist io___264 = { 0, 9, 0, 0, 0 };
                                static cilist io___265 = { 0, 10, 0, fmt_905, 0 };
                                static cilist io___266 = { 0, 12, 0, fmt_905, 0 };
                                static cilist io___267 = { 0, 10, 0, fmt_930, 0 };
                                static cilist io___268 = { 0, 9, 0, 0, 0 };
                                static cilist io___269 = { 0, 10, 0, fmt_905, 0 };
                                static cilist io___270 = { 0, 12, 0, fmt_905, 0 };
                                static cilist io___274 = { 0, 10, 0, fmt_891, 0 };
                                static cilist io___279 = { 0, 10, 0, fmt_1180, 0 };
                                static cilist io___280 = { 0, 10, 0, fmt_1200, 0 };
                                static cilist io___281 = { 0, 10, 0, fmt_1230, 0 };
                                static cilist io___282 = { 0, 10, 0, fmt_1240, 0 };
                                static cilist io___283 = { 0, 10, 0, fmt_1250, 0 };
                                static cilist io___284 = { 0, 10, 0, fmt_1260, 0 };
                                static cilist io___285 = { 0, 10, 0, fmt_1270, 0 };
                                static cilist io___286 = { 0, 10, 0, fmt_1310, 0 };
                                static cilist io___287 = { 0, 10, 0, fmt_850, 0 };
                                static cilist io___291 = { 0, 10, 0, fmt_1320, 0 };
                                static cilist io___292 = { 0, 10, 0, fmt_850, 0 };
                                static cilist io___293 = { 0, 10, 0, fmt_1350, 0 };
                                static cilist io___294 = { 0, 10, 0, fmt_1230, 0 };
                                static cilist io___295 = { 0, 10, 0, fmt_1240, 0 };
                                static cilist io___296 = { 0, 10, 0, fmt_1250, 0 };
                                static cilist io___297 = { 0, 10, 0, fmt_1260, 0 };
                                static cilist io___298 = { 0, 10, 0, fmt_1270, 0 };
                                static cilist io___299 = { 0, 10, 0, fmt_1310, 0 };
                                static cilist io___300 = { 0, 10, 0, fmt_850, 0 };
                                static cilist io___301 = { 0, 10, 0, fmt_1320, 0 };
                                static cilist io___302 = { 0, 10, 0, fmt_850, 0 };
                                static cilist io___303 = { 0, 12, 0, fmt_1435, 0 };
                                static cilist io___304 = { 0, 12, 0, fmt_1435, 0 };
                                static cilist io___305 = { 0, 12, 0, fmt_1435, 0 };
                                static cilist io___306 = { 0, 12, 0, fmt_1435, 0 };
                                static cilist io___307 = { 0, 12, 0, fmt_1435, 0 };
                                static cilist io___308 = { 0, 12, 0, fmt_1435, 0 };
                                static cilist io___309 = { 0, 6, 0, fmt_32, 0 };
                                static cilist io___310 = { 0, 6, 0, fmt_34, 0 };
                                static cilist io___311 = { 0, 6, 0, fmt_36, 0 };
                                static double result[100];


                                /*     SHIELDOSE-2, VERSION 2.10, 28 APR 94. */

                                /*        S.M. SELTZER */
                                /*        NATIONAL INSTITUTE OF STANDARDS AND TECHNOLOGY */
                                /*        GAITHERSBURG, MD 20899 */
                                /*        (301) 975-5552 */

                                /*        IDET = 1, AL DETECTOR */
                                /*               2, GRAPHITE DETECTOR */
                                /*               3, SI DETECTOR */
                                /*               4, AIR DETECTOR */
                                /*               5, BONE DETECTOR */
                                /*               6, CALCIUM FLUORIDE DETECTOR */
                                /*               7, GALLIUM ARSENIDE DETECTOR */
                                /*               8, LITHIUM FLUORIDE DETECTOR */
                                /*               9, SILICON DIOXIDE DETECTOR */
                                /*              10, TISSUE DETECTOR */
                                /*              11, WATER DETECTOR */

                                /*        INUC = 1, NO NUCLEAR ATTENUATION FOR PROTONS IN AL */
                                /*               2, NUCLEAR ATTENUATION, LOCAL CHARGED-SECONDARY ENERGY */
                                /*                     DEPOSITION */
                                /*               3, NUCLEAR ATTENUATION, LOCAL CHARGED-SECONDARY ENERGY */
                                /*                     DEPOSITION, AND APPROX EXPONENTIAL DISTRIBUTION OF */
                                /*                     NEUTRON DOSE */

                                /*        INCIDENT OMNIDIRECTIONAL FLUX IN /ENERGY/CM2/UNIT TIME */
                                /*             (SOLAR-FLARE FLUX IN /ENERGY/CM2). */

                                /*        EUNIT IS CONVERSION FACTOR FROM /ENERGY TO /MEV, */
                                /*             E.G., EUNIT = 1000 IF FLUX IS /KEV. */

                                /*        DURATN IS MISSION DURATION IN MULTIPLES OF UNIT TIME. */

                                /*     IMPLICIT DOUBLE PRECISION (A-H,O-Z) */
                                /*      CALL LOGO (VERSION) */

                                s_wsfe(&io___7);
                                //   e_wsfe();
                                //   s_rsfe(&io___8);
                                //   do_fio(&c__1, filenm, (ftnlen)40);
                                strcpy_s(filenm, 12, "example.inp");
                                //    e_rsfe();
                                o__1.oerr = 0;
                                o__1.ounit = 9;
                                o__1.ofnmlen = 40;
                                o__1.ofnm = filenm;
                                o__1.orl = 0;
                                o__1.osta = 0;
                                o__1.oacc = 0;
                                o__1.ofm = 0;
                                o__1.oblnk = 0;
                                strcpy_s(filenm, 12, "example.inp");
                                f_open(&o__1);
                                s_rsfe(&io___10);
                                do_fio(&c__1, prtfil, (ftnlen)40);
                                e_rsfe();
                                o__1.oerr = 0;
                                o__1.ounit = 10;
                                o__1.ofnmlen = 40;
                                o__1.ofnm = prtfil;
                                o__1.orl = 0;
                                o__1.osta = 0;
                                o__1.oacc = 0;
                                o__1.ofm = 0;
                                o__1.oblnk = 0;
                                f_open(&o__1);
                                s_rsfe(&io___12);
                                do_fio(&c__1, arrfil, (ftnlen)40);
                                e_rsfe();
                                o__1.oerr = 0;
                                o__1.ounit = 12;
                                o__1.ofnmlen = 40;
                                o__1.ofnm = arrfil;
                                o__1.orl = 0;
                                o__1.osta = 0;
                                o__1.oacc = 0;
                                o__1.ofm = 0;
                                o__1.oblnk = 0;
                                f_open(&o__1);
                                s_wsfe(&io___14);
                                do_fio(&c__1, version, (ftnlen)4);
                                e_wsfe();
                                s_wsfe(&io___15);
                                do_fio(&c__1, version, (ftnlen)4);
                                e_wsfe();
                                s_wsfe(&io___16);
                                do_fio(&c__1, filenm, (ftnlen)40);
                                e_wsfe();
                                s_wsfe(&io___17);
                                do_fio(&c__1, filenm, (ftnlen)40);
                                e_wsfe();
                                s_wsfe(&io___18);
                                do_fio(&c__1, prtfil, (ftnlen)40);
                                e_wsfe();
                                s_wsfe(&io___19);
                                do_fio(&c__1, prtfil, (ftnlen)40);
                                e_wsfe();
                                s_wsfe(&io___20);
                                do_fio(&c__1, arrfil, (ftnlen)40);
                                e_wsfe();
                                s_wsfe(&io___21);
                                do_fio(&c__1, arrfil, (ftnlen)40);
                                e_wsfe();
                                s_wsfe(&io___22);
                                e_wsfe();
                                s_rsle(&io___23);
                                do_lio(&c__3, &c__1, (char*)&idet, (ftnlen)sizeof(integer));
                                do_lio(&c__3, &c__1, (char*)&inuc, (ftnlen)sizeof(integer));
                                do_lio(&c__3, &c__1, (char*)&imax, (ftnlen)sizeof(integer));
                                do_lio(&c__3, &c__1, (char*)&iunt, (ftnlen)sizeof(integer));
                                e_rsle();
                                s_wsfe(&io___28);
                                do_fio(&c__1, (char*)&idet, (ftnlen)sizeof(integer));
                                do_fio(&c__1, (char*)&inuc, (ftnlen)sizeof(integer));
                                do_fio(&c__1, (char*)&imax, (ftnlen)sizeof(integer));
                                do_fio(&c__1, (char*)&iunt, (ftnlen)sizeof(integer));
                                e_wsfe();
                                s_wsfe(&io___29);
                                do_fio(&c__1, (char*)&idet, (ftnlen)sizeof(integer));
                                do_fio(&c__1, (char*)&imax, (ftnlen)sizeof(integer));
                                do_fio(&c__1, (char*)&inuc, (ftnlen)sizeof(integer));
                                e_wsfe();
                                inatt = 2;
                                if (inuc == 1) {
                                    inatt = 1;
                                }
                                inewt = 0;
                                if (inuc == 3) {
                                    inewt = 1;
                                }
                                s_wsle(&io___32);
                                //    do_lio(&c__9, &c__1, " Reading database and preparing base arrays......."
                                //	    "......", (ftnlen)56);
                                //    e_wsle();
                                //    s_wsle(&io___33);
                                 //   do_lio(&c__9, &c__1, "    Protons......................................."
                                //	    "......", (ftnlen)56);
                                e_wsle();
                                o__1.oerr = 0;
                                o__1.ounit = 11;
                                o__1.ofnmlen = 12;
                                o__1.ofnm = "PROTBAS2.DAT";
                                o__1.orl = 0;
                                o__1.osta = 0;
                                o__1.oacc = 0;
                                o__1.ofm = 0;
                                o__1.oblnk = 0;
                                f_open(&o__1);
                                s_rsfe(&io___34);
                                do_fio(&c__1, tag, (ftnlen)72);
                                e_rsfe();
                                s_rsle(&io___36);
                                do_lio(&c__3, &c__1, (char*)&mmaxp, (ftnlen)sizeof(integer));
                                do_lio(&c__3, &c__1, (char*)&kmaxp, (ftnlen)sizeof(integer));
                                do_lio(&c__3, &c__1, (char*)&nmaxp, (ftnlen)sizeof(integer));
                                do_lio(&c__3, &c__1, (char*)&lmaxp, (ftnlen)sizeof(integer));
                                do_lio(&c__3, &c__1, (char*)&imix, (ftnlen)sizeof(integer));
                                e_rsle();
                                s_rsle(&io___42);
                                i__1 = mmaxp;
                                for (m = 1; m <= i__1; ++m) {
                                    do_lio(&c__4, &c__1, (char*)&ep[m - 1], (ftnlen)sizeof(real));
                                }
                                e_rsle();
                                s_rsle(&io___45);
                                i__1 = mmaxp;
                                for (m = 1; m <= i__1; ++m) {
                                    do_lio(&c__4, &c__1, (char*)&rp[m - 1], (ftnlen)sizeof(real));
                                }
                                e_rsle();
                                s_rsle(&io___47);
                                i__1 = kmaxp;
                                for (k = 1; k <= i__1; ++k) {
                                    do_lio(&c__4, &c__1, (char*)&tepn[k - 1], (ftnlen)sizeof(real));
                                }
                                e_rsle();
                                s_rsle(&io___50);
                                i__1 = kmaxp;
                                for (k = 1; k <= i__1; ++k) {
                                    do_lio(&c__4, &c__1, (char*)&fepn[k - 1], (ftnlen)sizeof(real));
                                }
                                e_rsle();
                                s_rsle(&io___52);
                                i__1 = nmaxp;
                                for (n = 1; n <= i__1; ++n) {
                                    do_lio(&c__4, &c__1, (char*)&tp[n - 1], (ftnlen)sizeof(real));
                                }
                                e_rsle();
                                s_rsle(&io___55);
                                i__1 = lmaxp;
                                for (l = 1; l <= i__1; ++l) {
                                    do_lio(&c__4, &c__1, (char*)&zrp[l - 1], (ftnlen)sizeof(real));
                                }
                                e_rsle();
                                i__1 = nmaxp;
                                for (n = 1; n <= i__1; ++n) {
                                    for (i__ = 1; i__ <= 2; ++i__) {
                                        s_rsle(&io___59);
                                        i__2 = lmaxp;
                                        for (l = 1; l <= i__2; ++l) {
                                            do_lio(&c__4, &c__1, (char*)&dum[l - 1], (ftnlen)sizeof(real)
                                            );
                                        }
                                        e_rsle();
                                        if (i__ != inatt) {
                                            goto L38;
                                        }
                                        i__2 = lmaxp;
                                        for (l = 1; l <= i__2; ++l) {
                                            /* L390: */
                                            dalp[n + l * 49 - 50] = dum[l - 1];
                                        }
                                        dalp[n + lmaxp * 49 - 50] = dalp[n + (lmaxp - 1) * 49 - 50] /
                                            dalp[n + (lmaxp - 2) * 49 - 50] * dalp[n + (lmaxp - 1) *
                                            49 - 50];
                                    L38:
                                        ;
                                    }
                                    i__2 = imix;
                                    for (i__ = 1; i__ <= i__2; ++i__) {
                                        s_rsle(&io___62);
                                        i__3 = lmaxp;
                                        for (l = 1; l <= i__3; ++l) {
                                            do_lio(&c__4, &c__1, (char*)&dum[l - 1], (ftnlen)sizeof(real)
                                            );
                                        }
                                        e_rsle();
                                        if (i__ != idet) {
                                            goto L40;
                                        }
                                        i__3 = lmaxp;
                                        for (l = 1; l <= i__3; ++l) {
                                            /* L39: */
                                            dratp[n + l * 49 - 50] = dum[l - 1];
                                        }
                                    L40:
                                        ;
                                    }
                                    /* L50: */
                                }
                                cl__1.cerr = 0;
                                cl__1.cunit = 11;
                                cl__1.csta = 0;
                                f_clos(&cl__1);
                                i__1 = mmaxp;
                                for (m = 1; m <= i__1; ++m) {
                                    ep[m - 1] = log(ep[m - 1]);
                                    /* L60: */
                                    rp[m - 1] = log(rp[m - 1]);
                                }
                                i__1 = kmaxp;
                                for (k = 1; k <= i__1; ++k) {
                                    /* L65: */
                                    tepn[k - 1] = log(tepn[k - 1]);
                                }
                                i__1 = nmaxp;
                                for (n = 1; n <= i__1; ++n) {
                                    /* L70: */
                                    tp[n - 1] = log(tp[n - 1]);
                                }
                                scof_(&mmaxp, ep, rp, rpb, rpc, rpd);
                                scof_(&kmaxp, tepn, fepn, fepnb, fepnc, fepnd);
                                i__1 = lmaxp;
                                for (l = 1; l <= i__1; ++l) {
                                    i__2 = nmaxp;
                                    for (n = 1; n <= i__2; ++n) {
                                        /* L72: */
                                        dalp[n + l * 49 - 50] = log(dalp[n + l * 49 - 50]);
                                    }
                                    scof_(&nmaxp, tp, &dalp[l * 49 - 49], &dalpb[l * 49 - 49], &dalpc[l *
                                        49 - 49], &dalpd[l * 49 - 49]);
                                    /* L75: */
                                }
                                s_wsle(&io___73);
                                //   do_lio(&c__9, &c__1, "    Electrons and bremsstrahlung.................."
                               //	    "......", (ftnlen)56);
                               //    e_wsle();
                                o__1.oerr = 0;
                                o__1.ounit = 11;
                                o__1.ofnmlen = 12;
                                o__1.ofnm = "ELBRBAS2.DAT";
                                o__1.orl = 0;
                                o__1.osta = 0;
                                o__1.oacc = 0;
                                o__1.ofm = 0;
                                o__1.oblnk = 0;
                                f_open(&o__1);
                                s_rsfe(&io___74);
                                do_fio(&c__1, tag, (ftnlen)72);
                                e_rsfe();
                                s_rsle(&io___75);
                                do_lio(&c__3, &c__1, (char*)&mmaxe, (ftnlen)sizeof(integer));
                                do_lio(&c__3, &c__1, (char*)&nmaxe, (ftnlen)sizeof(integer));
                                do_lio(&c__3, &c__1, (char*)&lmaxs, (ftnlen)sizeof(integer));
                                do_lio(&c__3, &c__1, (char*)&lmaxe, (ftnlen)sizeof(integer));
                                do_lio(&c__3, &c__1, (char*)&lmaxt, (ftnlen)sizeof(integer));
                                do_lio(&c__3, &c__1, (char*)&lmaxb, (ftnlen)sizeof(integer));
                                do_lio(&c__3, &c__1, (char*)&imix, (ftnlen)sizeof(integer));
                                e_rsle();
                                nlene = nmaxe - nbege + 1;
                                s_rsle(&io___83);
                                i__1 = mmaxe;
                                for (m = 1; m <= i__1; ++m) {
                                    do_lio(&c__4, &c__1, (char*)&ee[m - 1], (ftnlen)sizeof(real));
                                }
                                e_rsle();
                                s_rsle(&io___85);
                                i__1 = mmaxe;
                                for (m = 1; m <= i__1; ++m) {
                                    do_lio(&c__4, &c__1, (char*)&re[m - 1], (ftnlen)sizeof(real));
                                }
                                e_rsle();
                                s_rsle(&io___87);
                                i__1 = mmaxe;
                                for (m = 1; m <= i__1; ++m) {
                                    do_lio(&c__4, &c__1, (char*)&ye[m - 1], (ftnlen)sizeof(real));
                                }
                                e_rsle();
                                s_rsle(&io___89);
                                i__1 = nmaxe;
                                for (n = 1; n <= i__1; ++n) {
                                    do_lio(&c__4, &c__1, (char*)&te[n - 1], (ftnlen)sizeof(real));
                                }
                                e_rsle();
                                s_rsle(&io___91);
                                i__1 = nmaxe;
                                for (n = 1; n <= i__1; ++n) {
                                    do_lio(&c__4, &c__1, (char*)&ar[n - 1], (ftnlen)sizeof(real));
                                }
                                e_rsle();
                                s_rsle(&io___93);
                                i__1 = nmaxe;
                                for (n = 1; n <= i__1; ++n) {
                                    do_lio(&c__4, &c__1, (char*)&rs[n - 1], (ftnlen)sizeof(real));
                                }
                                e_rsle();
                                s_rsle(&io___95);
                                i__1 = lmaxs;
                                for (l = 1; l <= i__1; ++l) {
                                    do_lio(&c__4, &c__1, (char*)&bs[l - 1], (ftnlen)sizeof(real));
                                }
                                e_rsle();
                                bs[lmaxs] = 2.f;
                                s_rsle(&io___97);
                                i__1 = lmaxe;
                                for (l = 1; l <= i__1; ++l) {
                                    do_lio(&c__4, &c__1, (char*)&zre[l - 1], (ftnlen)sizeof(real));
                                }
                                e_rsle();
                                s_rsle(&io___99);
                                i__1 = lmaxt;
                                for (l = 1; l <= i__1; ++l) {
                                    do_lio(&c__4, &c__1, (char*)&zs[l - 1], (ftnlen)sizeof(real));
                                }
                                e_rsle();
                                s_rsle(&io___101);
                                i__1 = lmaxb;
                                for (l = 1; l <= i__1; ++l) {
                                    do_lio(&c__4, &c__1, (char*)&zb[l - 1], (ftnlen)sizeof(real));
                                }
                                e_rsle();
                                i__1 = nmaxe;
                                for (n = 1; n <= i__1; ++n) {
                                    s_rsle(&io___103);
                                    i__2 = lmaxs;
                                    for (l = 1; l <= i__2; ++l) {
                                        do_lio(&c__4, &c__1, (char*)&dale[n + l * 14 - 15], (ftnlen)
                                            sizeof(real));
                                    }
                                    e_rsle();
                                    dale[n + (lmaxs + 1) * 14 - 15] = 1e-7f;
                                    s_rsle(&io___105);
                                    i__2 = lmaxt;
                                    for (l = 1; l <= i__2; ++l) {
                                        do_lio(&c__4, &c__1, (char*)&dalb[n + l * 14 - 15], (ftnlen)
                                            sizeof(real));
                                    }
                                    e_rsle();
                                    i__2 = imix;
                                    for (i__ = 1; i__ <= i__2; ++i__) {
                                        for (m = 1; m <= 2; ++m) {
                                            s_rsle(&io___107);
                                            i__3 = lmaxe;
                                            for (l = 1; l <= i__3; ++l) {
                                                do_lio(&c__4, &c__1, (char*)&dum[l - 1], (ftnlen)sizeof(
                                                    real));
                                            }
                                            e_rsle();
                                            if (i__ != idet) {
                                                goto L77;
                                            }
                                            i__3 = lmaxe;
                                            for (l = 1; l <= i__3; ++l) {
                                                /* L76: */
                                                drate[n + (l + m * 51) * 14 - 729] = dum[l - 1];
                                            }
                                        L77:
                                            s_rsle(&io___109);
                                            i__3 = lmaxb;
                                            for (l = 1; l <= i__3; ++l) {
                                                do_lio(&c__4, &c__1, (char*)&dum[l - 1], (ftnlen)sizeof(
                                                    real));
                                            }
                                            e_rsle();
                                            if (i__ != idet) {
                                                goto L80;
                                            }
                                            i__3 = lmaxb;
                                            for (l = 1; l <= i__3; ++l) {
                                                /* L78: */
                                                dratb[n + (l + m * 47) * 14 - 673] = dum[l - 1];
                                            }
                                        L80:
                                            ;
                                        }
                                        /* L90: */
                                    }
                                    /* L100: */
                                }
                                ++lmaxs;
                                cl__1.cerr = 0;
                                cl__1.cunit = 11;
                                cl__1.csta = 0;
                                f_clos(&cl__1);
                                i__1 = mmaxe;
                                for (m = 1; m <= i__1; ++m) {
                                    ee[m - 1] = log(ee[m - 1]);
                                    re[m - 1] = log(re[m - 1]);
                                    /* L110: */
                                    ye[m - 1] = log(ye[m - 1]);
                                }
                                i__1 = nmaxe;
                                for (n = 1; n <= i__1; ++n) {
                                    te[n - 1] = log(te[n - 1]);
                                    ar[n - 1] = log(ar[n - 1]);
                                    /* L120: */
                                    rs[n - 1] = log(rs[n - 1]);
                                }
                                i__1 = lmaxb;
                                for (l = 1; l <= i__1; ++l) {
                                    /* L130: */
                                    zb[l - 1] = log(zb[l - 1]);
                                }
                                scof_(&mmaxe, ee, re, reb, rec, red);
                                scof_(&mmaxe, ee, ye, yeb, yec, yed);
                                scof_(&nmaxe, te, ar, arb, arc, ard);
                                scof_(&nmaxe, te, rs, rsb, rsc, rsd);
                                i__1 = lmaxs;
                                for (l = 1; l <= i__1; ++l) {
                                    i__2 = nmaxe;
                                    for (n = nbege; n <= i__2; ++n) {
                                        /* L140: */
                                        dale[n + l * 14 - 15] = log(dale[n + l * 14 - 15]);
                                    }
                                    lcof_(&nlene, &te[nbege - 1], &dale[nbege + l * 14 - 15], &daleb[
                                        nbege + l * 14 - 15], &dalec[nbege + l * 14 - 15], &daled[
                                            nbege + l * 14 - 15]);
                                    /*     CALL SCOF (NLENE,TE(NBEGE),DALE(NBEGE,L),DALEB(NBEGE,L), */
                                    /*    1   DALEC(NBEGE,L),DALED(NBEGE,L)) */
                                    /* L150: */
                                }
                                i__1 = lmaxt;
                                for (l = 1; l <= i__1; ++l) {
                                    zs[l - 1] = log(zs[l - 1]);
                                    i__2 = nmaxe;
                                    for (n = nbege; n <= i__2; ++n) {
                                        /* L160: */
                                        dalb[n + l * 14 - 15] = log(dalb[n + l * 14 - 15]);
                                    }
                                    lcof_(&nlene, &te[nbege - 1], &dalb[nbege + l * 14 - 15], &dalbb[
                                        nbege + l * 14 - 15], &dalbc[nbege + l * 14 - 15], &dalbd[
                                            nbege + l * 14 - 15]);
                                    /*     CALL SCOF (NLENE,TE(NBEGE),DALB(NBEGE,L),DALBB(NBEGE,L), */
                                    /*    1   DALBC(NBEGE,L),DALBD(NBEGE,L)) */
                                    /* L170: */
                                }
                                /*     PRINT *,' Preparing base arrays for selected detector material...' */
                                i__1 = lmaxp;
                                for (l = 1; l <= i__1; ++l) {
                                    scof_(&nmaxp, tp, &dratp[l * 49 - 49], &dratpb[l * 49 - 49], &dratpc[
                                        l * 49 - 49], &dratpd[l * 49 - 49]);
                                    /* L220: */
                                }
                                for (m = 1; m <= 2; ++m) {
                                    i__1 = lmaxe;
                                    for (l = 1; l <= i__1; ++l) {
                                        lcof_(&nlene, &te[nbege - 1], &drate[nbege + (l + m * 51) * 14 -
                                            729], &drateb[nbege + (l + m * 51) * 14 - 729], &dratec[
                                                nbege + (l + m * 51) * 14 - 729], &drated[nbege + (l + m *
                                                    51) * 14 - 729]);
                                        /*     CALL SCOF (NLENE,TE(NBEGE),DRATE(NBEGE,L,M), */
                                        /*    1   DRATEB(NBEGE,L,M),DRATEC(NBEGE,L,M),DRATED(NBEGE,L,M)) */
                                        /* L230: */
                                    }
                                    /* L240: */
                                }
                                for (m = 1; m <= 2; ++m) {
                                    i__1 = lmaxb;
                                    for (l = 1; l <= i__1; ++l) {
                                        lcof_(&nlene, &te[nbege - 1], &dratb[nbege + (l + m * 47) * 14 -
                                            673], &dratbb[nbege + (l + m * 47) * 14 - 673], &dratbc[
                                                nbege + (l + m * 47) * 14 - 673], &dratbd[nbege + (l + m *
                                                    47) * 14 - 673]);
                                        /*     CALL SCOF (NLENE,TE(NBEGE),DRATB(NBEGE,L,M), */
                                        /*    1   DRATBB(NBEGE,L,M),DRATBC(NBEGE,L,M),DRATBD(NBEGE,L,M)) */
                                        /* L250: */
                                    }
                                    /* L260: */
                                }
                                switch (iunt) {
                                case 1:  goto L440;
                                case 2:  goto L470;
                                case 3:  goto L500;
                                }
                            L440:
                                s_wsfe(&io___138);
                                e_wsfe();
                                s_rsle(&io___139);
                                i__1 = imax;
                                for (i__ = 1; i__ <= i__1; ++i__) {
                                    do_lio(&c__4, &c__1, (char*)&zm[i__ - 1], (ftnlen)sizeof(real));
                                }
                                e_rsle();
                                s_wsfe(&io___141);
                                i__1 = imax;
                                for (i__ = 1; i__ <= i__1; ++i__) {
                                    do_fio(&c__1, (char*)&zm[i__ - 1], (ftnlen)sizeof(real));
                                }
                                e_wsfe();
                                i__1 = imax;
                                for (i__ = 1; i__ <= i__1; ++i__) {
                                    if (zm[i__ - 1] <= zmin / .0068580134999999993f) {
                                        zm[i__ - 1] = zmin / .0068580134999999993f;
                                    }
                                    z__[i__ - 1] = zm[i__ - 1] * .0068580134999999993f;
                                    /* L460: */
                                    zmm[i__ - 1] = z__[i__ - 1] * 3.7037037037037033f;
                                }
                                goto L530;
                            L470:
                                s_wsfe(&io___144);
                                e_wsfe();
                                s_rsle(&io___145);
                                i__1 = imax;
                                for (i__ = 1; i__ <= i__1; ++i__) {
                                    do_lio(&c__4, &c__1, (char*)&z__[i__ - 1], (ftnlen)sizeof(real));
                                }
                                e_rsle();
                                s_wsfe(&io___146);
                                i__1 = imax;
                                for (i__ = 1; i__ <= i__1; ++i__) {
                                    do_fio(&c__1, (char*)&z__[i__ - 1], (ftnlen)sizeof(real));
                                }
                                e_wsfe();
                                i__1 = imax;
                                for (i__ = 1; i__ <= i__1; ++i__) {
                                    if (z__[i__ - 1] <= zmin) {
                                        z__[i__ - 1] = zmin;
                                    }
                                    zm[i__ - 1] = z__[i__ - 1] / .0068580134999999993f;
                                    /* L490: */
                                    zmm[i__ - 1] = z__[i__ - 1] * 3.7037037037037033f;
                                }
                                goto L530;
                            L500:
                                s_wsfe(&io___147);
                                e_wsfe();
                                s_rsle(&io___148);
                                i__1 = imax;
                                for (i__ = 1; i__ <= i__1; ++i__) {
                                    do_lio(&c__4, &c__1, (char*)&zmm[i__ - 1], (ftnlen)sizeof(real));
                                }
                                e_rsle();
                                s_wsfe(&io___149);
                                i__1 = imax;
                                for (i__ = 1; i__ <= i__1; ++i__) {
                                    do_fio(&c__1, (char*)&zmm[i__ - 1], (ftnlen)sizeof(real));
                                }
                                e_wsfe();
                                i__1 = imax;
                                for (i__ = 1; i__ <= i__1; ++i__) {
                                    if (zmm[i__ - 1] <= zmin * 3.7037037037037033f) {
                                        zmm[i__ - 1] = zmin * 3.7037037037037033f;
                                    }
                                    z__[i__ - 1] = zmm[i__ - 1] / 3.7037037037037033f;
                                    /* L520: */
                                    zm[i__ - 1] = z__[i__ - 1] / .0068580134999999993f;
                                }
                            L530:
                                i__1 = imax;
                                for (i__ = 1; i__ <= i__1; ++i__) {
                                    /* L540: */
                                    zl[i__ - 1] = log(z__[i__ - 1]);
                                }
                                s_wsfe(&io___151);
                                i__1 = imax;
                                for (i__ = 1; i__ <= i__1; ++i__) {
                                    do_fio(&c__1, (char*)&z__[i__ - 1], (ftnlen)sizeof(real));
                                }
                                e_wsfe();
                                s_wsfe(&io___152);
                                e_wsfe();
                                s_rsle(&io___153);
                                do_lio(&c__4, &c__1, (char*)&emins, (ftnlen)sizeof(real));
                                do_lio(&c__4, &c__1, (char*)&emaxs, (ftnlen)sizeof(real));
                                do_lio(&c__4, &c__1, (char*)&eminp, (ftnlen)sizeof(real));
                                do_lio(&c__4, &c__1, (char*)&emaxp, (ftnlen)sizeof(real));
                                do_lio(&c__3, &c__1, (char*)&nptsp, (ftnlen)sizeof(integer));
                                do_lio(&c__4, &c__1, (char*)&emine, (ftnlen)sizeof(real));
                                do_lio(&c__4, &c__1, (char*)&emaxe, (ftnlen)sizeof(real));
                                do_lio(&c__3, &c__1, (char*)&nptse, (ftnlen)sizeof(integer));
                                e_rsle();
                                s_wsfe(&io___162);
                                do_fio(&c__1, (char*)&emins, (ftnlen)sizeof(real));
                                do_fio(&c__1, (char*)&emaxs, (ftnlen)sizeof(real));
                                do_fio(&c__1, (char*)&eminp, (ftnlen)sizeof(real));
                                do_fio(&c__1, (char*)&emaxp, (ftnlen)sizeof(real));
                                do_fio(&c__1, (char*)&nptsp, (ftnlen)sizeof(integer));
                                do_fio(&c__1, (char*)&emine, (ftnlen)sizeof(real));
                                do_fio(&c__1, (char*)&emaxe, (ftnlen)sizeof(real));
                                do_fio(&c__1, (char*)&nptse, (ftnlen)sizeof(integer));
                                e_wsfe();
                                eminu = dmin(eminp, emins);
                                emaxu = dmax(emaxp, emaxs);
                                dep = log(emaxu / eminu) / (real)(nptsp - 1);
                                eminul = log(eminu);
                                delp = dep / 3.f;
                                eindex_(&eminu, &dep, &nptsp, &emins, &emaxs, &nfstsb, &nlstsb, &nlensb);
                                eindex_(&eminu, &dep, &nptsp, &eminp, &emaxp, &nfstpb, &nlstpb, &nlenpb);
                                i__1 = nptsp;
                                for (np = 1; np <= i__1; ++np) {
                                    tpl[np - 1] = eminul + (real)(np - 1) * dep;
                                    tpp[np - 1] = exp(tpl[np - 1]);
                                    /* L570: */
                                }
                                s_wsfe(&io___177);
                                do_fio(&c__1, (char*)&tpp[nfstsb - 1], (ftnlen)sizeof(real));
                                do_fio(&c__1, (char*)&tpp[nlstsb - 1], (ftnlen)sizeof(real));
                                do_fio(&c__1, (char*)&tpp[nfstpb - 1], (ftnlen)sizeof(real));
                                do_fio(&c__1, (char*)&tpp[nlstpb - 1], (ftnlen)sizeof(real));
                                do_fio(&c__1, (char*)&nptsp, (ftnlen)sizeof(integer));
                                do_fio(&c__1, (char*)&emine, (ftnlen)sizeof(real));
                                do_fio(&c__1, (char*)&emaxe, (ftnlen)sizeof(real));
                                do_fio(&c__1, (char*)&nptse, (ftnlen)sizeof(integer));
                                e_wsfe();
                                s_wsfe(&io___178);
                                do_fio(&c__1, (char*)&tpp[nfstsb - 1], (ftnlen)sizeof(real));
                                do_fio(&c__1, (char*)&tpp[nlstsb - 1], (ftnlen)sizeof(real));
                                do_fio(&c__1, (char*)&tpp[nfstpb - 1], (ftnlen)sizeof(real));
                                do_fio(&c__1, (char*)&tpp[nlstpb - 1], (ftnlen)sizeof(real));
                                do_fio(&c__1, (char*)&nptsp, (ftnlen)sizeof(integer));
                                do_fio(&c__1, (char*)&emine, (ftnlen)sizeof(real));
                                do_fio(&c__1, (char*)&emaxe, (ftnlen)sizeof(real));
                                do_fio(&c__1, (char*)&nptse, (ftnlen)sizeof(integer));
                                e_wsfe();
                                s_wsle(&io___179);
                                //    do_lio(&c__9, &c__1, " Preparing mesh arrays to be integrated over spect"
                                //	    "ra....", (ftnlen)56);
                                //    e_wsle();
                                s_wsle(&io___180);
                                //    do_lio(&c__9, &c__1, "    Protons......................................."
                                //	    "......", (ftnlen)56);
                                //    e_wsle();
                                i__1 = nptsp;
                                for (np = 1; np <= i__1; ++np) {
                                    bspol_(&tpl[np - 1], &mmaxp, ep, rp, rpb, rpc, rpd, &ans);
                                    rinp = exp(ans);
                                    i__2 = lmaxp;
                                    for (l = 1; l <= i__2; ++l) {
                                        if (tpl[np - 1] < tp[nmaxp - 1]) {
                                            goto L590;
                                        }
                                        ans = dalp[nmaxp + l * 49 - 50];
                                        ansr = dratp[nmaxp + l * 49 - 50];
                                        goto L605;
                                    L590:
                                        if (tpl[np - 1] > tp[0]) {
                                            goto L600;
                                        }
                                        ans = dalp[l * 49 - 49];
                                        ansr = dratp[l * 49 - 49];
                                        goto L605;
                                    L600:
                                        bspol_(&tpl[np - 1], &nmaxp, tp, &dalp[l * 49 - 49], &dalpb[l *
                                            49 - 49], &dalpc[l * 49 - 49], &dalpd[l * 49 - 49], &ans);
                                        ansr = 1.f;
                                        if (idet == 1) {
                                            goto L605;
                                        }
                                        bspol_(&tpl[np - 1], &nmaxp, tp, &dratp[l * 49 - 49], &dratpb[l *
                                            49 - 49], &dratpc[l * 49 - 49], &dratpd[l * 49 - 49], &
                                            ansr);
                                    L605:
                                        din[l - 1] = ans + log(ansr);
                                        /* L610: */
                                    }
                                    enewt[np - 1] = 0.f;
                                    benmu = 0.f;
                                    if (inatt == 1) {
                                        goto L620;
                                    }
                                    if (tpl[np - 1] <= tepn[0]) {
                                        goto L615;
                                    }
                                    bspol_(&tpl[np - 1], &kmaxp, tepn, fepn, fepnb, fepnc, fepnd, &ans);
                                    enewt[np - 1] = tpp[np - 1] * ans;
                                L615:
                                    benmu = enewt[np - 1] * enmu;
                                L620:
                                    scof_(&lmaxp, zrp, din, dinb, dinc, dind);
                                    i__2 = imax;
                                    for (i__ = 1; i__ <= i__2; ++i__) {
                                        zrin = z__[i__ - 1] / rinp;
                                        if (zrin < zrp[lmaxp - 1]) {
                                            goto L640;
                                        }
                                        gp[np + i__ * 1001 - 1002] = 0.f;
                                        goto L645;
                                    L640:
                                        bspol_(&zrin, &lmaxp, zrp, din, dinb, dinc, dind, &ans);
                                        ans = exp(ans);
                                        gp[np + i__ * 1001 - 1002] = tpp[np - 1] * ans / rinp;
                                    L645:
                                        if (inewt == 1 && tpl[np - 1] > tepn[0]) {
                                            gp[np + i__ * 1001 - 1002] += benmu * exp(-enmu * z__[i__ - 1]
                                            );
                                        }
                                        /* L650: */
                                    }
                                    /* L660: */
                                }
                                s_wsle(&io___192);
                                //   do_lio(&c__9, &c__1, "    Electrons and bremsstrahlung.................."
                               //	    "......", (ftnlen)56);
                                //   e_wsle();
                                eminel = log(emine);
                                dee = (log(emaxe) - eminel) / (real)(nptse - 1);
                                dele = dee / 3.f;
                                i__1 = nptse;
                                for (ne = 1; ne <= i__1; ++ne) {
                                    tel[ne - 1] = eminel + (real)(ne - 1) * dee;
                                    tee[ne - 1] = exp(tel[ne - 1]);
                                    bspol_(&tel[ne - 1], &mmaxe, ee, re, reb, rec, red, &ans);
                                    rine[ne - 1] = exp(ans);
                                    bspol_(&tel[ne - 1], &nmaxe, te, rs, rsb, rsc, rsd, &ans);
                                    rins[ne - 1] = rine[ne - 1] * exp(ans);
                                    bspol_(&tel[ne - 1], &nmaxe, te, ar, arb, arc, ard, &ans);
                                    ares[ne - 1] = exp(ans);
                                    bspol_(&tel[ne - 1], &mmaxe, ee, ye, yeb, yec, yed, &ans);
                                    /* L670: */
                                    ylde[ne - 1] = exp(ans);
                                }
                                for (m = 1; m <= 2; ++m) {
                                    i__1 = nptse;
                                    for (ne = 1; ne <= i__1; ++ne) {
                                        i__2 = lmaxs;
                                        for (l = 1; l <= i__2; ++l) {
                                            if (tel[ne - 1] < te[nmaxe - 1]) {
                                                goto L680;
                                            }
                                            din[l - 1] = dale[nmaxe + l * 14 - 15];
                                            goto L700;
                                        L680:
                                            if (tel[ne - 1] > te[nbege - 1]) {
                                                goto L690;
                                            }
                                            din[l - 1] = dale[nbege + l * 14 - 15];
                                            goto L700;
                                        L690:
                                            bspol_(&tel[ne - 1], &nlene, &te[nbege - 1], &dale[nbege + l *
                                                14 - 15], &daleb[nbege + l * 14 - 15], &dalec[nbege
                                                + l * 14 - 15], &daled[nbege + l * 14 - 15], &din[l -
                                                1]);
                                        L700:
                                            ;
                                        }
                                        i__2 = lmaxe;
                                        for (l = 1; l <= i__2; ++l) {
                                            drin[l - 1] = 1.f;
                                            if (idet == 1 && m == 1) {
                                                goto L715;
                                            }
                                            if (tel[ne - 1] < te[nmaxe - 1]) {
                                                goto L710;
                                            }
                                            drin[l - 1] = drate[nmaxe + (l + m * 51) * 14 - 729];
                                            goto L715;
                                        L710:
                                            if (tel[ne - 1] > te[nbege - 1]) {
                                                goto L712;
                                            }
                                            drin[l - 1] = drate[nbege + (l + m * 51) * 14 - 729];
                                            goto L715;
                                        L712:
                                            bspol_(&tel[ne - 1], &nlene, &te[nbege - 1], &drate[nbege + (
                                                l + m * 51) * 14 - 729], &drateb[nbege + (l + m * 51)
                                                * 14 - 729], &dratec[nbege + (l + m * 51) * 14 - 729],
                                                &drated[nbege + (l + m * 51) * 14 - 729], &drin[l -
                                                1]);
                                            if (drin[l - 1] < 0.f) {
                                                drin[l - 1] = 0.f;
                                            }
                                        L715:
                                            ;
                                        }
                                        /*     CALL LCOF (LMAXS,BS,DIN,DINB,DINC,DIND) */
                                        scof_(&lmaxs, bs, din, dinb, dinc, dind);
                                        i__2 = imax;
                                        for (i__ = 1; i__ <= i__2; ++i__) {
                                            zrin = z__[i__ - 1] / rins[ne - 1];
                                            if (zrin < bs[lmaxs - 1]) {
                                                goto L730;
                                            }
                                            /* L720: */
                                            ge[ne + (i__ + m * 71) * 1001 - 72073] = 0.f;
                                            goto L740;
                                        L730:
                                            bspol_(&zrin, &lmaxs, bs, din, dinb, dinc, dind, &ans);
                                            ans = exp(ans);
                                            ge[ne + (i__ + m * 71) * 1001 - 72073] = tee[ne - 1] * ans *
                                                ares[ne - 1] / rins[ne - 1];
                                        L740:
                                            ;
                                        }
                                        /*     CALL LCOF (LMAXE,ZRE,DRIN,DINB,DINC,DIND) */
                                        scof_(&lmaxe, zre, drin, dinb, dinc, dind);
                                        i__2 = imax;
                                        for (i__ = 1; i__ <= i__2; ++i__) {
                                            zrin = z__[i__ - 1] / rine[ne - 1];
                                            if (zrin < zre[lmaxe - 1]) {
                                                goto L742;
                                            }
                                            ge[ne + (i__ + m * 71) * 1001 - 72073] *= drin[lmaxe - 1];
                                            goto L745;
                                        L742:
                                            bspol_(&zrin, &lmaxe, zre, drin, dinb, dinc, dind, &ansr);
                                            if (ansr < 0.f) {
                                                ansr = 0.f;
                                            }
                                            ge[ne + (i__ + m * 71) * 1001 - 72073] *= ansr;
                                        L745:
                                            ;
                                        }
                                        i__2 = lmaxt;
                                        for (l = 1; l <= i__2; ++l) {
                                            if (tel[ne - 1] < te[nmaxe - 1]) {
                                                goto L760;
                                            }
                                            din[l - 1] = dalb[nmaxe + l * 14 - 15];
                                            goto L780;
                                        L760:
                                            if (tel[ne - 1] > te[nbege - 1]) {
                                                goto L770;
                                            }
                                            din[l - 1] = dalb[nbege + l * 14 - 15];
                                            goto L780;
                                        L770:
                                            bspol_(&tel[ne - 1], &nlene, &te[nbege - 1], &dalb[nbege + l *
                                                14 - 15], &dalbb[nbege + l * 14 - 15], &dalbc[nbege
                                                + l * 14 - 15], &dalbd[nbege + l * 14 - 15], &din[l -
                                                1]);
                                        L780:
                                            ;
                                        }
                                        i__2 = lmaxb;
                                        for (l = 1; l <= i__2; ++l) {
                                            drin[l - 1] = 1.f;
                                            if (idet == 1 && m == 1) {
                                                goto L795;
                                            }
                                            if (tel[ne - 1] < te[nmaxe - 1]) {
                                                goto L790;
                                            }
                                            drin[l - 1] = dratb[nmaxe + (l + m * 47) * 14 - 673];
                                            goto L795;
                                        L790:
                                            if (tel[ne - 1] > te[nbege - 1]) {
                                                goto L792;
                                            }
                                            drin[l - 1] = dratb[nbege + (l + m * 47) * 14 - 673];
                                            goto L795;
                                        L792:
                                            bspol_(&tel[ne - 1], &nlene, &te[nbege - 1], &dratb[nbege + (
                                                l + m * 47) * 14 - 673], &dratbb[nbege + (l + m * 47)
                                                * 14 - 673], &dratbc[nbege + (l + m * 47) * 14 - 673],
                                                &dratbd[nbege + (l + m * 47) * 14 - 673], &drin[l -
                                                1]);
                                            if (drin[l - 1] < 0.f) {
                                                drin[l - 1] = 0.f;
                                            }
                                        L795:
                                            ;
                                        }
                                        lcof_(&lmaxt, zs, din, dinb, dinc, dind);
                                        /*     CALL SCOF (LMAXT,ZS,DIN,DINB,DINC,DIND) */
                                        i__2 = imax;
                                        for (i__ = 1; i__ <= i__2; ++i__) {
                                            zrinl = log(z__[i__ - 1] / rine[ne - 1]);
                                            bspol_(&zrinl, &lmaxt, zs, din, dinb, dinc, dind, &ans);
                                            ans = exp(ans);
                                            gb[ne + (i__ + m * 71) * 1001 - 72073] = tee[ne - 1] * ans *
                                                ylde[ne - 1] / rine[ne - 1];
                                            /* L800: */
                                        }
                                        lcof_(&lmaxb, zb, drin, dinb, dinc, dind);
                                        /*     CALL SCOF (LMAXB,ZB,DRIN,DINB,DINC,DIND) */
                                        i__2 = imax;
                                        for (i__ = 1; i__ <= i__2; ++i__) {
                                            if (zl[i__ - 1] < zb[lmaxb - 1]) {
                                                goto L810;
                                            }
                                            gb[ne + (i__ + m * 71) * 1001 - 72073] *= drin[lmaxb - 1];
                                            goto L812;
                                        L810:
                                            bspol_(&zl[i__ - 1], &lmaxb, zb, drin, dinb, dinc, dind, &
                                                ansr);
                                            if (ansr < 0.f) {
                                                ansr = 0.f;
                                            }
                                            gb[ne + (i__ + m * 71) * 1001 - 72073] *= ansr;
                                        L812:
                                            ;
                                        }
                                        /* L815: */
                                    }
                                    /* L820: */
                                }
                                s_wsle(&io___207);
                                //    do_lio(&c__9, &c__1, " Performing calculations for input spectra........"
                                //	    "......", (ftnlen)56);
                                 //   e_wsle();
                            L830:
                                s_wsfe(&io___208);
                                e_wsfe();
                                s_wsfe(&io___209);
                                e_wsfe();
                                i__1 = s_rsfe(&io___210);
                                if (i__1 != 0) {
                                    goto L1440;
                                }
                                i__1 = do_fio(&c__1, tag, (ftnlen)72);
                                if (i__1 != 0) {
                                    goto L1440;
                                }
                                i__1 = e_rsfe();
                                if (i__1 != 0) {
                                    goto L1440;
                                }
                                //s_wsfe(&io___211);
                                //do_fio(&c__1, tag, (ftnlen)72);
                                //e_wsfe();
                                s_wsfe(&io___212);
                                do_fio(&c__1, tag, (ftnlen)72);
                                e_wsfe();
                                s_wsfe(&io___213);
                                do_fio(&c__1, tag, (ftnlen)72);
                                e_wsfe();
                                s_wsfe(&io___214);
                                e_wsfe();
                                s_rsle(&io___215);
                                do_lio(&c__3, &c__1, (char*)&jsmax, (ftnlen)sizeof(integer));
                                do_lio(&c__3, &c__1, (char*)&jpmax, (ftnlen)sizeof(integer));
                                do_lio(&c__3, &c__1, (char*)&jemax, (ftnlen)sizeof(integer));
                                do_lio(&c__4, &c__1, (char*)&eunit, (ftnlen)sizeof(real));
                                do_lio(&c__4, &c__1, (char*)&duratn, (ftnlen)sizeof(real));
                                e_rsle();
                                s_wsfe(&io___221);
                                do_fio(&c__1, (char*)&jsmax, (ftnlen)sizeof(integer));
                                do_fio(&c__1, (char*)&jpmax, (ftnlen)sizeof(integer));
                                do_fio(&c__1, (char*)&jemax, (ftnlen)sizeof(integer));
                                do_fio(&c__1, (char*)&eunit, (ftnlen)sizeof(real));
                                do_fio(&c__1, (char*)&duratn, (ftnlen)sizeof(real));
                                e_wsfe();
                                s_wsfe(&io___222);
                                do_fio(&c__1, (char*)&jsmax, (ftnlen)sizeof(integer));
                                do_fio(&c__1, (char*)&jpmax, (ftnlen)sizeof(integer));
                                do_fio(&c__1, (char*)&jemax, (ftnlen)sizeof(integer));
                                do_fio(&c__1, (char*)&eunit, (ftnlen)sizeof(real));
                                do_fio(&c__1, (char*)&duratn, (ftnlen)sizeof(real));
                                e_wsfe();
                                if (duratn <= 0.f) {
                                    duratn = 1.f;
                                }
                                deltas = radcon * delp / 4.f;
                                deltap = duratn * radcon * delp / 4.f;
                                deltae = duratn * radcon * dele / 4.f;
                                if (eunit <= 0.f) {
                                    eunit = 1.f;
                                }
                                isol = 2;
                                if (jsmax < 3) {
                                    goto L900;
                                }
                                isol = 1;
                                s_wsfe(&io___227);
                                e_wsfe();
                                s_rsle(&io___228);
                                i__1 = jsmax;
                                for (j = 1; j <= i__1; ++j) {
                                    do_lio(&c__4, &c__1, (char*)&eps[j - 1], (ftnlen)sizeof(real));
                                }
                                e_rsle();
                                s_wsfe(&io___231);
                                i__1 = jsmax;
                                for (j = 1; j <= i__1; ++j) {
                                    do_fio(&c__1, (char*)&eps[j - 1], (ftnlen)sizeof(real));
                                }
                                e_wsfe();
                                s_wsfe(&io___232);
                                i__1 = jsmax;
                                for (j = 1; j <= i__1; ++j) {
                                    do_fio(&c__1, (char*)&eps[j - 1], (ftnlen)sizeof(real));
                                }
                                e_wsfe();
                                s_wsfe(&io___233);
                                e_wsfe();
                                s_rsle(&io___234);
                                i__1 = jsmax;
                                for (j = 1; j <= i__1; ++j) {
                                    do_lio(&c__4, &c__1, (char*)&s[j - 1], (ftnlen)sizeof(real));
                                }
                                e_rsle();
                                s_wsfe(&io___236);
                                i__1 = jsmax;
                                for (j = 1; j <= i__1; ++j) {
                                    do_fio(&c__1, (char*)&s[j - 1], (ftnlen)sizeof(real));
                                }
                                e_wsfe();
                                s_wsfe(&io___237);
                                i__1 = jsmax;
                                for (j = 1; j <= i__1; ++j) {
                                    do_fio(&c__1, (char*)&s[j - 1], (ftnlen)sizeof(real));
                                }
                                e_wsfe();
                                nlens = nlensb;
                                nfsts = nfstsb;
                                nlsts = nlstsb;
                                spectr_(&jsmax, eps, s, &eunit, &eminu, &dep, &nptsp, &nfsts, &nlsts, &
                                    nlens, tpp, tpl, sol);
                                s_wsfe(&io___242);
                                do_fio(&c__1, (char*)&tpp[nfsts - 1], (ftnlen)sizeof(real));
                                do_fio(&c__1, (char*)&tpp[nlsts - 1], (ftnlen)sizeof(real));
                                do_fio(&c__1, (char*)&nlens, (ftnlen)sizeof(integer));
                                e_wsfe();
                                i__1 = nlsts;
                                for (np = nfsts; np <= i__1; ++np) {
                                    /* L892: */
                                    g[np - 1] = sol[np - 1] * enewt[np - 1];
                                }
                                integ_(&delp, &g[nfsts - 1], &nlens, &eneut);
                                i__1 = nlsts;
                                for (np = nfsts; np <= i__1; ++np) {
                                    /* L894: */
                                    g[np - 1] = sol[np - 1] * tpp[np - 1];
                                }
                                integ_(&delp, &g[nfsts - 1], &nlens, &eav);
                                eneut /= eav;
                                s_wsfe(&io___246);
                                do_fio(&c__1, (char*)&eneut, (ftnlen)sizeof(real));
                                e_wsfe();
                            L900:
                                itrp = 2;
                                if (jpmax < 3) {
                                    goto L920;
                                }
                                itrp = 1;
                                s_wsfe(&io___248);
                                e_wsfe();
                                s_rsle(&io___249);
                                i__1 = jpmax;
                                for (j = 1; j <= i__1; ++j) {
                                    do_lio(&c__4, &c__1, (char*)&eps[j - 1], (ftnlen)sizeof(real));
                                }
                                e_rsle();
                                s_wsfe(&io___250);
                                i__1 = jpmax;
                                for (j = 1; j <= i__1; ++j) {
                                    do_fio(&c__1, (char*)&eps[j - 1], (ftnlen)sizeof(real));
                                }
                                e_wsfe();
                                s_wsfe(&io___251);
                                i__1 = jpmax;
                                for (j = 1; j <= i__1; ++j) {
                                    do_fio(&c__1, (char*)&eps[j - 1], (ftnlen)sizeof(real));
                                }
                                e_wsfe();
                                s_wsfe(&io___252);
                                e_wsfe();
                                s_rsle(&io___253);
                                i__1 = jpmax;
                                for (j = 1; j <= i__1; ++j) {
                                    do_lio(&c__4, &c__1, (char*)&s[j - 1], (ftnlen)sizeof(real));
                                }
                                e_rsle();
                                s_wsfe(&io___254);
                                i__1 = jpmax;
                                for (j = 1; j <= i__1; ++j) {
                                    do_fio(&c__1, (char*)&s[j - 1], (ftnlen)sizeof(real));
                                }
                                e_wsfe();
                                s_wsfe(&io___255);
                                i__1 = jpmax;
                                for (j = 1; j <= i__1; ++j) {
                                    do_fio(&c__1, (char*)&s[j - 1], (ftnlen)sizeof(real));
                                }
                                e_wsfe();
                                nlenp = nlenpb;
                                nfstp = nfstpb;
                                nlstp = nlstpb;
                                spectr_(&jpmax, eps, s, &eunit, &eminu, &dep, &nptsp, &nfstp, &nlstp, &
                                    nlenp, tpp, tpl, spg);
                                s_wsfe(&io___260);
                                do_fio(&c__1, (char*)&tpp[nfstp - 1], (ftnlen)sizeof(real));
                                do_fio(&c__1, (char*)&tpp[nlstp - 1], (ftnlen)sizeof(real));
                                do_fio(&c__1, (char*)&nlenp, (ftnlen)sizeof(integer));
                                e_wsfe();
                                i__1 = nlstp;
                                for (np = nfstp; np <= i__1; ++np) {
                                    /* L912: */
                                    g[np - 1] = spg[np - 1] * enewt[np - 1];
                                }
                                integ_(&delp, &g[nfstp - 1], &nlenp, &eneut);
                                i__1 = nlstp;
                                for (np = nfstp; np <= i__1; ++np) {
                                    /* L914: */
                                    g[np - 1] = spg[np - 1] * tpp[np - 1];
                                }
                                integ_(&delp, &g[nfstp - 1], &nlenp, &eav);
                                eneut /= eav;
                                s_wsfe(&io___261);
                                do_fio(&c__1, (char*)&eneut, (ftnlen)sizeof(real));
                                e_wsfe();
                            L920:
                                ilec = 2;
                                if (jemax < 3) {
                                    goto L940;
                                }
                                ilec = 1;
                                s_wsfe(&io___263);
                                e_wsfe();
                                s_rsle(&io___264);
                                i__1 = jemax;
                                for (j = 1; j <= i__1; ++j) {
                                    do_lio(&c__4, &c__1, (char*)&eps[j - 1], (ftnlen)sizeof(real));
                                }
                                e_rsle();
                                s_wsfe(&io___265);
                                i__1 = jemax;
                                for (j = 1; j <= i__1; ++j) {
                                    do_fio(&c__1, (char*)&eps[j - 1], (ftnlen)sizeof(real));
                                }
                                e_wsfe();
                                s_wsfe(&io___266);
                                i__1 = jemax;
                                for (j = 1; j <= i__1; ++j) {
                                    do_fio(&c__1, (char*)&eps[j - 1], (ftnlen)sizeof(real));
                                }
                                e_wsfe();
                                s_wsfe(&io___267);
                                e_wsfe();
                                s_rsle(&io___268);
                                i__1 = jemax;
                                for (j = 1; j <= i__1; ++j) {
                                    do_lio(&c__4, &c__1, (char*)&s[j - 1], (ftnlen)sizeof(real));
                                }
                                e_rsle();
                                s_wsfe(&io___269);
                                i__1 = jemax;
                                for (j = 1; j <= i__1; ++j) {
                                    do_fio(&c__1, (char*)&s[j - 1], (ftnlen)sizeof(real));
                                }
                                e_wsfe();
                                s_wsfe(&io___270);
                                i__1 = jemax;
                                for (j = 1; j <= i__1; ++j) {
                                    do_fio(&c__1, (char*)&s[j - 1], (ftnlen)sizeof(real));
                                }
                                e_wsfe();
                                nlene = nptse;
                                nfste = 1;
                                nlste = nptse;
                                spectr_(&jemax, eps, s, &eunit, &emine, &dee, &nptse, &nfste, &nlste, &
                                    nlene, tee, tel, seg);
                                s_wsfe(&io___274);
                                do_fio(&c__1, (char*)&tee[nfste - 1], (ftnlen)sizeof(real));
                                do_fio(&c__1, (char*)&tee[nlste - 1], (ftnlen)sizeof(real));
                                do_fio(&c__1, (char*)&nlene, (ftnlen)sizeof(integer));
                                e_wsfe();
                            L940:
                                switch (isol) {
                                case 1:  goto L980;
                                case 2:  goto L950;
                                }
                            L950:
                                i__1 = nlsts;
                                for (np = nfsts; np <= i__1; ++np) {
                                    /* L960: */
                                    sol[np - 1] = 0.f;
                                }
                                for (j = 1; j <= 2; ++j) {
                                    i__1 = imax;
                                    for (i__ = 1; i__ <= i__1; ++i__) {
                                        /* L970: */
                                        dosol[i__ + j * 71 - 72] = 0.f;
                                    }
                                }
                                goto L1010;
                            L980:
                                i__1 = imax;
                                for (i__ = 1; i__ <= i__1; ++i__) {
                                    i__2 = nlsts;
                                    for (np = nfsts; np <= i__2; ++np) {
                                        /* L990: */
                                        g[np - 1] = sol[np - 1] * gp[np + i__ * 1001 - 1002];
                                    }
                                    integ_(&deltas, &g[nfsts - 1], &nlens, &dosol[i__ - 1]);
                                    /* L1000: */
                                }
                                sphere_(zl, dosol, &imax, &dosol[71]);
                            L1010:
                                switch (itrp) {
                                case 1:  goto L1050;
                                case 2:  goto L1020;
                                }
                            L1020:
                                i__1 = nlstp;
                                for (np = nfstp; np <= i__1; ++np) {
                                    /* L1030: */
                                    spg[np - 1] = 0.f;
                                }
                                for (j = 1; j <= 2; ++j) {
                                    i__1 = imax;
                                    for (i__ = 1; i__ <= i__1; ++i__) {
                                        /* L1040: */
                                        dosp[i__ + j * 71 - 72] = 0.f;
                                    }
                                }
                                goto L1080;
                            L1050:
                                i__1 = imax;
                                for (i__ = 1; i__ <= i__1; ++i__) {
                                    i__2 = nlstp;
                                    for (np = nfstp; np <= i__2; ++np) {
                                        /* L1060: */
                                        g[np - 1] = spg[np - 1] * gp[np + i__ * 1001 - 1002];
                                    }
                                    integ_(&deltap, &g[nfstp - 1], &nlenp, &dosp[i__ - 1]);
                                    /* L1070: */
                                }
                                sphere_(zl, dosp, &imax, &dosp[71]);
                            L1080:
                                switch (ilec) {
                                case 1:  goto L1110;
                                case 2:  goto L1090;
                                }
                            L1090:
                                for (j = 1; j <= 2; ++j) {
                                    for (m = 1; m <= 2; ++m) {
                                        i__1 = imax;
                                        for (i__ = 1; i__ <= i__1; ++i__) {
                                            dose[i__ + (m + (j << 1)) * 71 - 214] = 0.f;
                                            /* L1100: */
                                            dosb[i__ + (m + (j << 1)) * 71 - 214] = 0.f;
                                        }
                                    }
                                }
                                goto L1160;
                            L1110:
                                for (m = 1; m <= 2; ++m) {
                                    i__1 = imax;
                                    for (i__ = 1; i__ <= i__1; ++i__) {
                                        i__2 = nlste;
                                        for (ne = nfste; ne <= i__2; ++ne) {
                                            g[ne - 1] = seg[ne - 1] * ge[ne + (i__ + m * 71) * 1001 -
                                                72073];
                                            /* L1120: */
                                            spg[ne - 1] = seg[ne - 1] * gb[ne + (i__ + m * 71) * 1001 -
                                                72073];
                                        }
                                        integ_(&deltae, &g[nfste - 1], &nlene, &dose[i__ + (m + 2) * 71 -
                                            214]);
                                        integ_(&deltae, &spg[nfste - 1], &nlene, &dosb[i__ + (m + 2) * 71
                                            - 214]);
                                        /* L1130: */
                                    }
                                    switch (m) {
                                    case 1:  goto L1140;
                                    case 2:  goto L1150;
                                    }
                                L1140:
                                    sphere_(zl, &dose[(m + 2) * 71 - 213], &imax, &dose[(m + 4) * 71 -
                                        213]);
                                    sphere_(zl, &dosb[(m + 2) * 71 - 213], &imax, &dosb[(m + 4) * 71 -
                                        213]);
                                L1150:
                                    ;
                                }
                            L1160:
                                j = 1;
                                for (m = 2; m >= 1; --m) {
                                    switch (m) {
                                    case 1:  goto L1190;
                                    case 2:  goto L1170;
                                    }
                                L1170:
                                    s_wsfe(&io___279);
                                    e_wsfe();
                                    goto L1210;
                                L1190:
                                    s_wsfe(&io___280);
                                    e_wsfe();
                                L1210:
                                    s_wsfe(&io___281);
                                    do_fio(&c__1, det + (idet - 1 << 3), (ftnlen)8);
                                    e_wsfe();
                                    if (inatt == 1) {
                                        s_wsfe(&io___282);
                                        e_wsfe();
                                    }
                                    if (inatt == 2) {
                                        s_wsfe(&io___283);
                                        e_wsfe();
                                    }
                                    if (inatt == 2 && inewt == 0) {
                                        s_wsfe(&io___284);
                                        e_wsfe();
                                    }
                                    if (inatt == 2 && inewt == 1) {
                                        s_wsfe(&io___285);
                                        e_wsfe();
                                    }
                                    s_wsfe(&io___286);
                                    e_wsfe();
                                    s_wsfe(&io___287);
                                    e_wsfe();
                                    i__1 = imax;

                                    for (i__ = 1; i__ <= i__1; ++i__) {
                                        doseb = dose[i__ + (m + (j << 1)) * 71 - 214] + dosb[i__ + (m + (
                                            j << 1)) * 71 - 214];
                                        dosebp = doseb + dosp[i__ + j * 71 - 72];
                                        dost = dosebp + dosol[i__ + j * 71 - 72];
                                        result[i__] = (double)dost;
                                        s_wsfe(&io___291);
                                        do_fio(&c__1, (char*)&zm[i__ - 1], (ftnlen)sizeof(real));
                                        do_fio(&c__1, (char*)&zmm[i__ - 1], (ftnlen)sizeof(real));
                                        do_fio(&c__1, (char*)&z__[i__ - 1], (ftnlen)sizeof(real));
                                        do_fio(&c__1, (char*)&dose[i__ + (m + (j << 1)) * 71 - 214], (
                                            ftnlen)sizeof(real));
                                        do_fio(&c__1, (char*)&dosb[i__ + (m + (j << 1)) * 71 - 214], (
                                            ftnlen)sizeof(real));
                                        do_fio(&c__1, (char*)&doseb, (ftnlen)sizeof(real));
                                        do_fio(&c__1, (char*)&dosp[i__ + j * 71 - 72], (ftnlen)sizeof(
                                            real));
                                        do_fio(&c__1, (char*)&dosol[i__ + j * 71 - 72], (ftnlen)sizeof(
                                            real));
                                        do_fio(&c__1, (char*)&dosebp, (ftnlen)sizeof(real));
                                        do_fio(&c__1, (char*)&dost, (ftnlen)sizeof(real));
                                        e_wsfe();
                                        if ((real)(i__ / 10) == (real)i__ * .1f) {
                                            s_wsfe(&io___292);
                                            e_wsfe();
                                        }
                                        /* L1330: */
                                    }
                                    /* L1340: */
                                }
                                j = 2;
                                m = 1;
                                s_wsfe(&io___293);
                                e_wsfe();
                                s_wsfe(&io___294);
                                do_fio(&c__1, det + (idet - 1 << 3), (ftnlen)8);
                                e_wsfe();
                                if (inatt == 1) {
                                    s_wsfe(&io___295);
                                    e_wsfe();
                                }
                                if (inatt == 2) {
                                    s_wsfe(&io___296);
                                    e_wsfe();
                                }
                                if (inatt == 2 && inewt == 0) {
                                    s_wsfe(&io___297);
                                    e_wsfe();
                                }
                                if (inatt == 2 && inewt == 1) {
                                    s_wsfe(&io___298);
                                    e_wsfe();
                                }
                                s_wsfe(&io___299);
                                e_wsfe();
                                s_wsfe(&io___300);
                                e_wsfe();
                                i__1 = imax;
                                for (i__ = 1; i__ <= i__1; ++i__) {
                                    doseb = dose[i__ + (m + (j << 1)) * 71 - 214] + dosb[i__ + (m + (j <<
                                        1)) * 71 - 214];
                                    dosebp = doseb + dosp[i__ + j * 71 - 72];
                                    dost = dosebp + dosol[i__ + j * 71 - 72];
                                    result[i__] = (double)dost;
                                    //printf("%f", result[i__]);
                                    s_wsfe(&io___301);
                                    do_fio(&c__1, (char*)&zm[i__ - 1], (ftnlen)sizeof(real));
                                    do_fio(&c__1, (char*)&zmm[i__ - 1], (ftnlen)sizeof(real));
                                    do_fio(&c__1, (char*)&z__[i__ - 1], (ftnlen)sizeof(real));
                                    do_fio(&c__1, (char*)&dose[i__ + (m + (j << 1)) * 71 - 214], (ftnlen)
                                        sizeof(real));
                                    do_fio(&c__1, (char*)&dosb[i__ + (m + (j << 1)) * 71 - 214], (ftnlen)
                                        sizeof(real));
                                    do_fio(&c__1, (char*)&doseb, (ftnlen)sizeof(real));
                                    do_fio(&c__1, (char*)&dosp[i__ + j * 71 - 72], (ftnlen)sizeof(real));
                                    do_fio(&c__1, (char*)&dosol[i__ + j * 71 - 72], (ftnlen)sizeof(real))
                                        ;
                                    do_fio(&c__1, (char*)&dosebp, (ftnlen)sizeof(real));
                                    do_fio(&c__1, (char*)&dost, (ftnlen)sizeof(real));
                                    e_wsfe();
                                    if ((real)(i__ / 10) == (real)i__ * .1f) {
                                        s_wsfe(&io___302);
                                        e_wsfe();
                                    }
                                    /* L1410: */
                                }
                                for (j = 1; j <= 2; ++j) {
                                    s_wsfe(&io___303);
                                    i__1 = imax;
                                    for (i__ = 1; i__ <= i__1; ++i__) {
                                        do_fio(&c__1, (char*)&dosol[i__ + j * 71 - 72], (ftnlen)sizeof(
                                            real));
                                    }
                                    e_wsfe();
                                    s_wsfe(&io___304);
                                    i__1 = imax;
                                    for (i__ = 1; i__ <= i__1; ++i__) {
                                        do_fio(&c__1, (char*)&dosp[i__ + j * 71 - 72], (ftnlen)sizeof(
                                            real));
                                    }
                                    e_wsfe();
                                    /* L1437: */
                                }
                                for (m = 2; m >= 1; --m) {
                                    s_wsfe(&io___305);
                                    i__1 = imax;
                                    for (i__ = 1; i__ <= i__1; ++i__) {
                                        do_fio(&c__1, (char*)&dose[i__ + (m + 2) * 71 - 214], (ftnlen)
                                            sizeof(real));
                                    }
                                    e_wsfe();
                                    s_wsfe(&io___306);
                                    i__1 = imax;
                                    for (i__ = 1; i__ <= i__1; ++i__) {
                                        do_fio(&c__1, (char*)&dosb[i__ + (m + 2) * 71 - 214], (ftnlen)
                                            sizeof(real));
                                    }
                                    e_wsfe();
                                    /* L1438: */
                                }
                                s_wsfe(&io___307);
                                i__1 = imax;
                                for (i__ = 1; i__ <= i__1; ++i__) {
                                    do_fio(&c__1, (char*)&dose[i__ + 141], (ftnlen)sizeof(real));
                                }
                                e_wsfe();
                                s_wsfe(&io___308);
                                i__1 = imax;
                                for (i__ = 1; i__ <= i__1; ++i__) {
                                    do_fio(&c__1, (char*)&dosb[i__ + 141], (ftnlen)sizeof(real));
                                }
                                e_wsfe();
                                goto L830;
                            L1440:
                                /*s_wsfe(&io___309);
                                do_fio(&c__1, filenm, (ftnlen)40);
                                e_wsfe();
                                s_wsfe(&io___310);
                                do_fio(&c__1, prtfil, (ftnlen)40);
                                e_wsfe();
                                s_wsfe(&io___311);
                                do_fio(&c__1, arrfil, (ftnlen)40);
                                e_wsfe();*/

                                /*      print 1500,CHAR(27) */
                                /* 1500 format (' ',A1,'[0m') */
                                //    s_stop("", (ftnlen)0);
                                return result[1];
} /* MAIN__ */

/*     SUBROUTINE EINDEX (EMINB,DE,NPTS,EMIN,EMAX,NFST,NLST,NLEN), 28 APR 94. */
/* Subroutine */ int eindex_(real* eminb, real* de, integer* npts, real* emin,
    real* emax, integer* nfst, integer* nlst, integer* nlen)
{
    /* Builtin functions */
    double log(doublereal);

    /*     IMPLICIT DOUBLE PRECISION (A-H,O-Z) */
    *nfst = log(*emin / *eminb) / *de + .5f;
    ++(*nfst);
    if (*nfst < 1) {
        *nfst = 1;
    }
    *nlst = log(*emax / *eminb) / *de + .5f;
    ++(*nlst);
    if (*nlst > * npts) {
        *nlst = *npts;
    }
    *nlen = *nlst - *nfst + 1;
    return 0;
} /* eindex_ */

/*     SUBROUTINE SPECTR (JMAX,EPS,S,EUNIT,EMINB,DEL,NPTS,NFST,NLST,NLEN, */
/*    1   T,TL,SP), 28 APR 94. */
/* Subroutine */ int spectr_(integer* jmax, real* eps, real* s, real* eunit,
    real* eminb, real* del, integer* npts, integer* nfst, integer* nlst,
    integer* nlen, real* t, real* tl, real* sp)
{
    /* Initialized data */

    static real emc2 = 938.27231f;
    static real arglim = 85.f;

    /* Format strings */
    static char fmt_60[] = "(/\002    INT SPEC    EAV(MeV)\002)";
    static char fmt_70[] = "(1pe12.4,0pf12.5)";

    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal), exp(doublereal), log(doublereal);
    integer s_wsfe(cilist*), e_wsfe(void), do_fio(integer*, char*, ftnlen);

    /* Local variables */
    static real g[1001];
    static integer j, n;
    static real p, arg, ans, sin__, bcof[301], ccof[301], dcof[301], beta,
        ebar;
    extern /* Subroutine */ int scof_(integer*, real*, real*, real*, real
        *, real*);
    static real alpha, delta;
    extern /* Subroutine */ int integ_(real*, real*, integer*, real*),
        bspol_(real*, integer*, real*, real*, real*, real*, real*,
            real*), eindex_(real*, real*, integer*, real*, real*,
                integer*, integer*, integer*);

    /* Fortran I/O blocks */
    static cilist io___328 = { 0, 10, 0, fmt_60, 0 };
    static cilist io___329 = { 0, 10, 0, fmt_70, 0 };


    /*     IMPLICIT DOUBLE PRECISION (A-H,O-Z) */
    /*           ARGLIM = 700.0 (DOUBLE PRECISION), 85.0 (SINGLE PRECISION) */
        /* Parameter adjustments */
    --sp;
    --tl;
    --t;
    --s;
    --eps;

    /* Function Body */
    delta = *del / 3.f;
    if (eps[1] > 0.f) {
        goto L20;
    }
    alpha = s[1];
    beta = s[2];
    if (beta <= 0.f) {
        beta = 1.f;
    }
    beta /= alpha;
    i__1 = *nlst;
    for (n = *nfst; n <= i__1; ++n) {
        sp[n] = 0.f;
        g[n - 1] = 0.f;
        if (s[3] <= 0.f) {
            goto L6;
        }
        p = sqrt(t[n] * (t[n] + emc2 * 2.f));
        arg = p / alpha;
        if (arg > arglim) {
            goto L10;
        }
        sp[n] = t[n] * beta * ((t[n] + emc2) / p) * exp(-arg);
        goto L8;
    L6:
        arg = t[n] / alpha;
        if (arg > arglim) {
            goto L10;
        }
        sp[n] = t[n] * beta * exp(-arg);
    L8:
        g[n - 1] = t[n] * sp[n];
    L10:
        ;
    }
    goto L50;
L20:
    eindex_(eminb, del, npts, &eps[1], &eps[*jmax], nfst, nlst, nlen);
    i__1 = *jmax;
    for (j = 1; j <= i__1; ++j) {
        eps[j] = log(eps[j]);
        /* L30: */
        s[j] = log(*eunit * s[j]);
    }
    scof_(jmax, &eps[1], &s[1], bcof, ccof, dcof);
    i__1 = *nlst;
    for (n = *nfst; n <= i__1; ++n) {
        bspol_(&tl[n], jmax, &eps[1], &s[1], bcof, ccof, dcof, &ans);
        sp[n] = t[n] * exp(ans);
        /* L40: */
        g[n - 1] = t[n] * sp[n];
    }
L50:
    integ_(&delta, &sp[*nfst], nlen, &sin__);
    integ_(&delta, &g[*nfst - 1], nlen, &ebar);
    ebar /= sin__;
    s_wsfe(&io___328);
    e_wsfe();
    s_wsfe(&io___329);
    do_fio(&c__1, (char*)&sin__, (ftnlen)sizeof(real));
    do_fio(&c__1, (char*)&ebar, (ftnlen)sizeof(real));
    e_wsfe();
    return 0;
} /* spectr_ */

/*     SUBROUTINE SPHERE (ZL,DOSE,IMAX,DOSPH), 26 JAN 93. */
/* Subroutine */ int sphere_(real* zl, real* dose, integer* imax, real* dosph)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static integer i__;
    static real bcof[71], ccof[71], dcof[71];
    extern /* Subroutine */ int scof_(integer*, real*, real*, real*, real
        *, real*);
    static real dosl[71];
    static integer imix, imix1;

    /*     IMPLICIT DOUBLE PRECISION (A-H,O-Z) */
        /* Parameter adjustments */
    --dosph;
    --dose;
    --zl;

    /* Function Body */
    i__1 = *imax;
    for (i__ = 1; i__ <= i__1; ++i__) {
        if (dose[i__] <= 0.f) {
            goto L20;
        }
        /* L10: */
        dosl[i__ - 1] = log(dose[i__]);
    }
    i__ = *imax + 1;
L20:
    imix = i__ - 1;
    if (imix < 3) {
        goto L40;
    }
    scof_(&imix, &zl[1], dosl, bcof, ccof, dcof);
    bcof[imix - 1] = bcof[imix - 2] + (ccof[imix - 2] * 2.f + dcof[imix - 2] *
        3.f * (zl[imix] - zl[imix - 1])) * (zl[imix] - zl[imix - 1]);
    i__1 = imix;
    for (i__ = 1; i__ <= i__1; ++i__) {
        /* L30: */
        dosph[i__] = dose[i__] * (1.f - bcof[i__ - 1]);
    }
L40:
    imix1 = imix + 1;
    if (imix1 > * imax) {
        return 0;
    }
    i__1 = *imax;
    for (i__ = imix1; i__ <= i__1; ++i__) {
        /* L50: */
        dosph[i__] = 0.f;
    }
    return 0;
} /* sphere_ */

/* Subroutine */ int scof_(integer* n, real* x, real* y, real* b, real* c__,
    real* d__)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer j;
    static real r__, s;
    static integer n1, jr;

    /*        REINSCH ALGORITHM, VIA MJB, 22 FEB 83 */
    /*        Y(S)=((D(J)*(X-X(J))+C(J))*(X-X(J))+B(J))*(X-X(J))+Y(J) */
    /*             FOR X BETWEEN X(J) AND X(J+1) */
    /*     IMPLICIT DOUBLE PRECISION (A-H,O-Z) */
        /* Parameter adjustments */
    --d__;
    --c__;
    --b;
    --y;
    --x;

    /* Function Body */
    n1 = *n - 1;
    s = 0.f;
    i__1 = n1;
    for (j = 1; j <= i__1; ++j) {
        d__[j] = x[j + 1] - x[j];
        r__ = (y[j + 1] - y[j]) / d__[j];
        c__[j] = r__ - s;
        /* L10: */
        s = r__;
    }
    s = 0.f;
    r__ = 0.f;
    c__[1] = 0.f;
    c__[*n] = 0.f;
    i__1 = n1;
    for (j = 2; j <= i__1; ++j) {
        c__[j] += r__ * c__[j - 1];
        b[j] = (x[j - 1] - x[j + 1]) * 2.f - r__ * s;
        s = d__[j];
        /* L20: */
        r__ = s / b[j];
    }
    for (jr = n1; jr >= 2; --jr) {
        /* L30: */
        c__[jr] = (d__[jr] * c__[jr + 1] - c__[jr]) / b[jr];
    }
    i__1 = n1;
    for (j = 1; j <= i__1; ++j) {
        s = d__[j];
        r__ = c__[j + 1] - c__[j];
        d__[j] = r__ / s;
        c__[j] *= 3.f;
        /* L40: */
        b[j] = (y[j + 1] - y[j]) / s - (c__[j] + r__) * s;
    }
    return 0;
} /* scof_ */

/* Subroutine */ int bspol_(real* s, integer* n, real* x, real* y, real* b,
    real* c__, real* d__, real* t)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static real q;
    static integer ml, mu, mlb, mub, mav, idir, idir1;

    /*        BINARY SEARCH, X ASCENDING OR DESCENDING */
    /*     IMPLICIT DOUBLE PRECISION (A-H,O-Z) */
        /* Parameter adjustments */
    --d__;
    --c__;
    --b;
    --y;
    --x;

    /* Function Body */
    if (x[1] > x[*n]) {
        goto L10;
    }
    idir = 0;
    mlb = 0;
    mub = *n;
    goto L20;
L10:
    idir = 1;
    mlb = *n;
    mub = 0;
L20:
    idir1 = idir - 1;
    if (*s >= x[mub + idir]) {
        goto L60;
    }
    if (*s <= x[mlb - idir1]) {
        goto L70;
    }
    ml = mlb;
    mu = mub;
    goto L40;
L30:
    if ((i__1 = mu - ml, abs(i__1)) <= 1) {
        goto L80;
    }
L40:
    mav = (ml + mu) / 2;
    if (*s < x[mav]) {
        goto L50;
    }
    ml = mav;
    goto L30;
L50:
    mu = mav;
    goto L30;
L60:
    mu = mub + idir + idir1;
    goto L90;
L70:
    mu = mlb - idir - idir1;
    goto L90;
L80:
    mu += idir1;
L90:
    q = *s - x[mu];
    *t = ((d__[mu] * q + c__[mu]) * q + b[mu]) * q + y[mu];
    return 0;
} /* bspol_ */

/* Subroutine */ int lcof_(integer* nmax, real* x, real* f, real* b, real*
    c__, real* d__)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer n;

    /*        26 JAN 93.  SIMPLE LINEAR INTERPOLATION */
    /*     IMPLICIT DOUBLE PRECISION (A-H,O-Z) */
        /* Parameter adjustments */
    --d__;
    --c__;
    --b;
    --f;
    --x;

    /* Function Body */
    i__1 = *nmax - 1;
    for (n = 1; n <= i__1; ++n) {
        b[n] = (f[n + 1] - f[n]) / (x[n + 1] - x[n]);
        c__[n] = 0.f;
        d__[n] = 0.f;
        /* L10: */
    }
    return 0;
} /* lcof_ */

/* Subroutine */ int integ_(real* delta, real* g, integer* n, real* result)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer k, nl1, nl2;
    static real sig6, sum2, sum4, sigma;

    /*          INCLUDES N=1 */
    /*     IMPLICIT DOUBLE PRECISION (A-H,O-Z) */
        /* Parameter adjustments */
    --g;

    /* Function Body */
    nl1 = *n - 1;
    nl2 = *n - 2;
    if ((real)(*n) - (real)(*n / 2) * 2.f <= 0.f) {
        goto L100;
    }
    else {
        goto L10;
    }
L10:
    if (*n - 1 <= 0) {
        goto L15;
    }
    else {
        goto L20;
    }
L15:
    sigma = 0.f;
    goto L70;
L20:
    if (*n - 3 <= 0) {
        goto L30;
    }
    else {
        goto L40;
    }
L30:
    sigma = g[1] + g[2] * 4.f + g[3];
    goto L70;
L40:
    sum4 = 0.f;
    i__1 = nl1;
    for (k = 2; k <= i__1; k += 2) {
        /* L50: */
        sum4 += g[k];
    }
    sum2 = 0.f;
    i__1 = nl2;
    for (k = 3; k <= i__1; k += 2) {
        /* L60: */
        sum2 += g[k];
    }
    sigma = g[1] + sum4 * 4.f + sum2 * 2.f + g[*n];
L70:
    *result = *delta * sigma;
    return 0;
L100:
    if (*n - 2 <= 0) {
        goto L110;
    }
    else {
        goto L120;
    }
L110:
    sigma = (g[1] + g[2]) * 1.5f;
    goto L70;
L120:
    if (*n - 4 <= 0) {
        goto L130;
    }
    else {
        goto L140;
    }
L130:
    sigma = (g[1] + g[2] * 3.f + g[3] * 3.f + g[4]) * 1.125f;
    goto L70;
L140:
    if (*n - 6 <= 0) {
        goto L150;
    }
    else {
        goto L160;
    }
L150:
    sigma = g[1] + g[2] * 3.875f + g[3] * 2.625f + g[4] * 2.625f + g[5] *
        3.875f + g[6];
    goto L70;
L160:
    if (*n - 8 <= 0) {
        goto L170;
    }
    else {
        goto L180;
    }
L170:
    sigma = g[1] + g[2] * 3.875f + g[3] * 2.625f + g[4] * 2.625f + g[5] *
        3.875f + g[6] * 2.f + g[7] * 4.f + g[8];
    goto L70;
L180:
    sig6 = g[1] + g[2] * 3.875f + g[3] * 2.625f + g[4] * 2.625f + g[5] *
        3.875f + g[6];
    sum4 = 0.f;
    i__1 = nl1;
    for (k = 7; k <= i__1; k += 2) {
        /* L190: */
        sum4 += g[k];
    }
    sum2 = 0.f;
    i__1 = nl2;
    for (k = 8; k <= i__1; k += 2) {
        /* L200: */
        sum2 += g[k];
    }
    sigma = sig6 + g[6] + sum4 * 4.f + sum2 * 2.f + g[*n];
    goto L70;
} /* integ_ */



/* Subroutine */ int cls_(void)
{
    /* Format strings */
    static char fmt_10[] = "(\002 \002,a1,\002[2J\002)";

    /* Builtin functions */
    integer s_wsfe(cilist*), do_fio(integer*, char*, ftnlen), e_wsfe(void);

    /* Fortran I/O blocks */
    static cilist io___378 = { 0, 6, 0, fmt_10, 0 };


    /*        8 SEP 88. */
    s_wsfe(&io___378);
    do_fio(&c__1, "\033", (ftnlen)1);
    e_wsfe();
    return 0;
} /* cls_ */