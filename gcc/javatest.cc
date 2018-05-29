#include <iostream>
#include "hrtree_headers.h"

using namespace std;

#ifdef GLOBAL
int main(){
    const char* blockfile = "rtree.bin";
    const char* datafile = "vector.data";
    
    Cache* cache = new Cache(10, 1024);
    RTree* tree = new RTree(blockfile, 1024, cache, 2);
    tree->bulkload(datafile, 32);
    printf("Vector count: %d.\r\n", tree->num_of_data );
    float point[] = {135, 35};
    SortedLinList* list = new SortedLinList();
    tree->NNQuery(point, list);
    printf("Vector count:.\r\n");

    printf("Vector count: %d.\r\n", list->get_num());

    delete tree;
    delete cache;
    return 0;
}
#endif

#ifdef HILBERT_NEXTA
int main(){
    long corner0[] = {123, 456};
    long corner1[] = {321, 654};
    long point0[] = {156, 546};
    bitmask_t point1[] = {210, 320};
    
    int stat =  hilbert_nextinbox(2, 8, 32, 1, corner0, corner1, point1);

    if (stat)
      for (int i = 0; i < 2; ++i)
        printf("%u ", corner0[i]);
    else
      printf("No such point");
    printf("\n");

    return 0;
}
#endif
#ifdef HILBERT_POINT_TEMP
int main(){
    double corner1[] = {789.789, 654.654};
    char* x = reinterpret_cast<char *>(corner1);
    propogateIEEEBits(0, 8, x , 1083, 0);
    cout << corner1[0] << "\t" << corner1[1] << endl;
    return 0;
}
#endif
#ifdef HILBERT_POINT_IEEE
int main(){
    int ndims=2, nbits=32;
    double corner0[] = {-123.123, 456.456};
    double corner1[] = {789.789, 654.654};
    double cornerlo[ndims];
    double cornerhi[ndims];
    double work[ndims];
    int k;

    for (int i = 0; i < ndims; i++) {
      cornerlo[i] = corner0[i];
      work[i] = corner1[i];
    }

    int result = hilbert_ieee_box_pt(ndims, 1, cornerlo, work);
    for (k = 0; k < ndims; ++k)
      printf("%20lf", cornerlo[k]);
    printf("\t %4u\n", result);

    for (int i = 0; i < ndims; i++) {
      work[i] = corner0[i];
      cornerhi[i] = corner1[i];
    }
    result = hilbert_ieee_box_pt(ndims, 0, work, cornerhi);
    for (k = 0; k < ndims; ++k)
      printf("%20lf", cornerhi[k]);
    printf("\t %4u\n", result);
    return 0;
}
#endif
#ifdef HILBERT_POINT
int main(){
    int ndims=2, nbits=32;
    long corner0[] = {123, 456};
    long corner1[] = {789, 654};
    long cornerlo[ndims];
    long cornerhi[ndims];
    long work[ndims];
    int k;

    for (int i = 0; i < ndims; i++) {
      cornerlo[i] = corner0[i];
      work[i] = corner1[i];
    }

    int result = hilbert_box_pt(ndims, 8, 32, 1, cornerlo, work);
    for (k = 0; k < ndims; ++k)
      printf("%8u", cornerlo[k]);
    printf("\t %4u\n", result);

    for (int i = 0; i < ndims; i++) {
      work[i] = corner0[i];
      cornerhi[i] = corner1[i];
    }
    result = hilbert_box_pt(ndims, 8, 32, 0, work, cornerhi);
    for (k = 0; k < ndims; ++k)
      printf("%8u", cornerhi[k]);
    printf("\t %4u\n", result);
    return 0;
}
#endif

#ifdef TEST_IEEE_VTX2
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#define maxDim (8*sizeof(bitmask_t))
typedef double key_tt;

unsigned g_nDims;

int cmp(void const* c1p, void const* c2p)
{
  double const* x = reinterpret_cast<double const*>(c1p);
  double const* y = reinterpret_cast<double const*>(c2p);
  return hilbert_ieee_cmp(g_nDims, x, y);
}

int main()
{
  key_tt corner0[maxDim], corner1[maxDim];
  key_tt cornerlo[maxDim], cornerhi[maxDim], work[maxDim];
  typedef key_tt array_t[maxDim];
  array_t* array;

  unsigned nDims, i, k;

  printf( "Enter nDims: " );
  scanf( "%d", &nDims);
  if ( nDims == 0 )
    return 0;

  //for (i = 0; i < 10000; ++i)
  for (i = 0; i < 1; ++i)
  {
    for (k = 0; k < nDims; ++k)
    {
      corner0[k] = 2.*drand48() - 1.;
      corner1[k] = 2.*drand48() - 1.;
    }
    
    /* start temp test */
    corner0[0] = -0.207070;
    corner0[1] = -0.293328;
    corner1[0] = 0.680971;
    corner1[1] = -0.106833;
    /* end temp test */

    /* find first corner */
    for (k = 0; k < nDims; ++k)
    {
      cornerlo[k] = corner0[k];
      work[k] = corner1[k];
    }

    int result1 =  hilbert_ieee_box_vtx(nDims, 1, cornerlo, work);

    /* find last corner */
    for (k = 0; k < nDims; ++k)
    {
      work[k] = corner0[k];
      cornerhi[k] = corner1[k];
    }

    int result2 = hilbert_ieee_box_vtx(nDims, 0, work, cornerhi);

    printf("Predicted corner0:  ");
    for (k = 0; k < nDims; ++k)
      printf("%20lf", corner0[k]);
    printf("\n");

    printf("Predicted corner1:  ");
    for (k = 0; k < nDims; ++k)
      printf("%20lf", corner1[k]);
    printf("\n");

    printf("Predicted cornerlo: ");
    for (k = 0; k < nDims; ++k)
      printf("%20lf", cornerlo[k]);
    printf("\t%d\n", result1);

    printf("Predicted cornerhi: ");
    for (k = 0; k < nDims; ++k)
      printf("%20lf", cornerhi[k]);
    printf("\t%d\n", result2);


    array = (array_t*) malloc(maxDim*sizeof(key_tt) << nDims);
    for (k = 0; k < (1<<nDims); ++k)
    {
      unsigned j;
      key_tt* eltk = &array[k][0];
      for (j = 0; j < nDims; ++j)
      {
        key_tt* src = ((k>>j)&1)? corner1: corner0;
        eltk[j] = src[j];
      }
    }

    g_nDims = nDims;
    qsort(array, (1<<nDims), maxDim*sizeof(key_tt), cmp);

    for (k = 0; k < (1<<nDims); k += (1 << nDims) - 1)
    {
      unsigned j;
      int mismatch = 0;
      key_tt* eltk = &array[k][0];
      for (j = 0; j < nDims & !mismatch; ++j)
      {
        mismatch = (eltk[j] != ((k==0)? cornerlo: cornerhi)[j]);
      }
      assert (!mismatch);
    }
    free((char*)array);
  }
  return 0;
}

#endif
#ifdef HILBERT_BOX
int main(){
    int ndims=3, nbits=21;
    long corner0[] = {436246, 88525, 3380};
    long corner1[] = {446859, 108961, 391371};
    long cornerlo[ndims];
    long cornerhi[ndims];
    long work[ndims];
    int k;

    for (int i = 0; i < ndims; i++) {
      cornerlo[i] = corner0[i];
      work[i] = corner1[i];
    }

    int result = hilbert_box_vtx(ndims, 8, 32, 1, cornerlo, work);
    for (k = 0; k < ndims; ++k)
      printf("%8u", cornerlo[k]);
    printf("\t %4u\n", result);

    for (int i = 0; i < ndims; i++) {
      work[i] = corner0[i];
      cornerhi[i] = corner1[i];
    }
    result = hilbert_box_vtx(ndims, 8, 32, 0, work, cornerhi);
    for (k = 0; k < ndims; ++k)
      printf("%8u", cornerhi[k]);
    printf("\t %4u\n", result);
    return 0;
}
#endif

#ifdef HILBERT_C2I
int main(){
    bitmask_t coord[2];
    coord[0] = 353145208L;
    coord[1] = 1161980015L;
    bitmask_t index = hilbert_c2i(2, 32, coord);
    cout << index << endl;
    return 0;
}
#endif

#ifdef HILBERT_I2C
int main(){
    bitmask_t index = 1234567890123456789L;
    bitmask_t coord[4];
    hilbert_i2c(2, 32, index, coord);
    return 0;
}
#endif

#ifdef TEST_READBIT
int main(){
/*
    bitmask_t coord1[] = {353145208L, 1161980015L};
    bitmask_t coord2[] = {353145308L, 1161982015L};
    int ans = hilbert_cmp(2, sizeof(coord1[0]), 32, coord1, coord2);
    cout << ans << endl;
    ans = hilbert_cmp(2, sizeof(coord1[0]), 32, coord2, coord1);
    cout << ans << endl;
    ans = hilbert_cmp(2, sizeof(coord1[0]), 32, coord2, coord2);
    cout << ans << endl;

    printf( "cmp r = 0x%0*x; r1 = 0x%0*x, ans = %2d\n", 67/4, 9547923, 67/4, 8934793, ans );
    printf( "cmp r = 0x%0*x; r1 = 0x%0*x\n", 67/4, 100, 67/4, 2000);
*/
        bitmask_t coord1[64];
        bitmask_t coord2[64];
        bitmask_t index = -9223372036854775808L;
        bitmask_t idx = 9223372036854774807L;
        index = -9223372036854775808L;
        idx = 9223372036854774807L;
        hilbert_i2c(64 , 1, index, coord1);
        hilbert_i2c(64 , 1, idx, coord2);
        int result = hilbert_cmp(64, 8, 1, coord1, coord2);
cout << "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4611686018427387904 -4611686018427387904 " << endl;
cout << " 0 0 4 12 28 28 28 28 28 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 540 4611686018427388444 4611686018427388444  " << endl;
//				cout << endl << "3,0,2 2,0,2 1,0,2 1,0,2 1,0,2 1,0,2 1,0,2 1,0,2 1,0,2 1,0,2 1,0,2 1,0,2 1,0,2 1,0,2 1,0,2 1,0,2 1,0,2 1,0,2 1,0,2 1,0,2 1,0,2 1,0,2 0,0,2 1,0,2 2,0,2 3,0,2 1,0,2 2,0,2 1,0,2 1,0,2 2,0,2 0,0,2" << endl;
				cout << result << endl;
    return 0;
}
#endif

#ifdef TEST_IEEEC

int main(){
  double c1[] = {-0.326012, -0.498183};
  double c2[] = {-0.455216, -0.258825};
  double c3[] = {-0.469326, -0.0559034};
  cout << hilbert_ieee_cmp(2, c1, c2) << endl;
  cout << hilbert_ieee_cmp(2, c3, c1) << endl;
  cout << hilbert_ieee_cmp(2, c1, c1) << endl;
//  double coord1[] = {2147483654.0, 2147483618.0};
//  double coord2[] = {2147483654.0, 2147483617.0};
//  cout << hilbert_ieee_cmp(2, coord1, coord2) << endl;
//  cout << hilbert_ieee_cmp(2, coord2, coord1) << endl;
//  cout << hilbert_ieee_cmp(2, coord2, coord2) << endl;
  //long c1[] = {1, 3, 1};
  //long c2[] = {1, 1, 1};
  //cout << hilbert_cmp(3, 8, 21, c1, c2) << endl;
//  double c1[] = {1.0, 0.0, 1.0};
//  double c2[] = {0.0, 0.0, 1.0};
//  cout << hilbert_ieee_cmp(3, c1, c2) << endl;
//  cout << hilbert_ieee_cmp(3, coord2, coord1) << endl;
//  cout << hilbert_ieee_cmp(3, coord1, coord1) << endl;
//    bitmask_t coord1[] = {353145208L, 1161980015L};
//    bitmask_t coord2[] = {353145308L, 1161982015L};
//    int ans = hilbert_cmp(2, sizeof(coord1[0]), 32, coord1, coord2);
//    bitmask_t coord1[] = {353145208L, 1161980015L};
//    bitmask_t coord2[] = {353145308L, 1161982015L};
//    int ans = hilbert_cmp(2, sizeof(coord1[0]), 32, coord1, coord2);
  return 0;
}
#endif

#ifdef IEEE_TEST
int cmp(const void* xv, const void* yv)
{
  double const* x = reinterpret_cast<double const*>(xv);
  double const* y = reinterpret_cast<double const*>(yv);
  /* return hilbert_cmp(2, 8, 64, x, y); */
  return hilbert_ieee_cmp(2, x, y);
}
int main() {
  double *coords;
  coords = (double*) malloc(2*100*sizeof(double));
  coords[0]=0.1165652789091054;
  coords[1]=-0.082186671764504;
  coords[2]=-0.27499427517730035;
  coords[3]=-0.2945363873602148;
  coords[4]=0.31115489741307023;
  coords[5]=0.13542625136995057;
  coords[6]=0.06357896888974757;
  coords[7]=0.03998065293921094;
  coords[8]=-0.4369681244164153;
  coords[9]=-0.07134933009278421;
  coords[10]=-0.3864652577212728;
  coords[11]=-0.4524794279094364;
  coords[12]=-0.27583230606542963;
  coords[13]=-0.1946671953858733;
  coords[14]=0.03243603563832809;
  coords[15]=0.21445495760259503;
  coords[16]=-0.17696199338141894;
  coords[17]=0.29554306959063836;
  coords[18]=0.26606432382029366;
  coords[19]=0.4867882453591481;
  coords[20]=0.3387974263595426;
  coords[21]=-0.4349771951936672;
  coords[22]=-0.1788299490206522;
  coords[23]=0.0878112442478498;
  coords[24]=-0.3144654668190491;
  coords[25]=-0.10777794427051968;
  coords[26]=-0.49124440762210175;
  coords[27]=-0.1953533203277017;
  coords[28]=-0.20603337143449885;
  coords[29]=-0.17768639992261082;
  coords[30]=0.3289416829545442;
  coords[31]=-0.405330159068347;
  coords[32]=0.1593792178590343;
  coords[33]=-0.1366740020917735;
  coords[34]=0.10303848659104031;
  coords[35]=-0.007330929996177349;
  coords[36]=0.3623103524498894;
  coords[37]=0.047546230273148415;
  coords[38]=0.04926334235648455;
  coords[39]=0.33574504156082685;
  coords[40]=-0.13757445542772873;
  coords[41]=-0.4543993128314613;
  coords[42]=-0.25438344055174644;
  coords[43]=-0.23467040600929345;
  coords[44]=-0.18377601687868395;
  coords[45]=-0.2088698460621159;
  coords[46]=-0.17325169951264752;
  coords[47]=0.3143380934968458;
  coords[48]=-0.2243375494766473;
  coords[49]=0.2507335581465092;
  coords[50]=0.32121831825031755;
  coords[51]=-0.35764781034428705;
  coords[52]=0.387478014448541;
  coords[53]=0.056861897325444444;
  coords[54]=-0.24862787119911478;
  coords[55]=0.1428621425215859;
  coords[56]=-0.36806311850985896;
  coords[57]=0.21060218888172932;
  coords[58]=-0.009873592822404964;
  coords[59]=0.21398172483203282;
  coords[60]=-0.2850690618321319;
  coords[61]=-0.043201021592477296;
  coords[62]=0.3972773761408649;
  coords[63]=0.45768001967841554;
  coords[64]=0.4587898933869905;
  coords[65]=0.23571871625223217;
  coords[66]=-0.057484137038071426;
  coords[67]=0.45245887680040286;
  coords[68]=-0.07512850711813612;
  coords[69]=-0.17652990995534112;
  coords[70]=-0.13684679232243857;
  coords[71]=-0.4694636904290226;
  coords[72]=0.04493792082006043;
  coords[73]=-0.11234866812249855;
  coords[74]=-0.0983673072100184;
  coords[75]=-0.31702031325347224;
  coords[76]=0.2385217167227044;
  coords[77]=0.20520314115848282;
  coords[78]=0.32122085420030544;
  coords[79]=0.3148649497122098;
  coords[80]=0.34993018143541066;
  coords[81]=0.21620399725470774;
  coords[82]=-0.2921558535005101;
  coords[83]=0.07742138322070802;
  coords[84]=-0.24876834486770383;
  coords[85]=0.42035202831210283;
  coords[86]=0.03185624569030132;
  coords[87]=0.4152535397356363;
  coords[88]=0.12651373992692438;
  coords[89]=0.36213577393836116;
  coords[90]=-0.43440840308082074;
  coords[91]=-0.21588440433029843;
  coords[92]=0.2778973835993124;
  coords[93]=-0.40379621130878374;
  coords[94]=0.030881296230048094;
  coords[95]=-0.28165444209259394;
  coords[96]=0.27523726755807354;
  coords[97]=-0.42379907011971674;
  coords[98]=-0.21472226431005548;
  coords[99]=0.05694994872569359;
  coords[100]=0.4091334277732872;
  coords[101]=0.3683469994422648;
  coords[102]=-0.05344543402703483;
  coords[103]=-0.3265380158739164;
  coords[104]=-0.2991751801011041;
  coords[105]=0.4552910910802188;
  coords[106]=0.10211891646165827;
  coords[107]=-0.002920178575925969;
  coords[108]=0.1727306820023653;
  coords[109]=-0.3337170075961765;
  coords[110]=-0.3467067633838299;
  coords[111]=0.46122325117718443;
  coords[112]=-0.3101491600044487;
  coords[113]=0.4464606816203738;
  coords[114]=-0.09600575895550545;
  coords[115]=-0.22946898383426118;
  coords[116]=-0.0550685920962497;
  coords[117]=0.26106084604509583;
  coords[118]=-0.28084298565870425;
  coords[119]=-0.4615469691199636;
  coords[120]=0.07924783633369858;
  coords[121]=-0.08038301171057127;
  coords[122]=-0.04967766783215333;
  coords[123]=0.4149748906276163;
  coords[124]=-0.4939477424487453;
  coords[125]=7.177609006049357E-4;
  coords[126]=-0.14587428959822235;
  coords[127]=0.2976449631584802;
  coords[128]=-0.2648223328584227;
  coords[129]=0.35965590660479796;
  coords[130]=0.4139973213704976;
  coords[131]=-0.40104305669370843;
  coords[132]=0.38498668682215964;
  coords[133]=-0.09069434634491624;
  coords[134]=-0.18379867061166122;
  coords[135]=-0.1707263515915275;
  coords[136]=0.11179787294345422;
  coords[137]=0.47882486090087406;
  coords[138]=0.32841498822454207;
  coords[139]=-0.053139392935889096;
  coords[140]=-0.11288190163144252;
  coords[141]=-0.4247144908875694;
  coords[142]=0.4823662141854248;
  coords[143]=0.02523110805165263;
  coords[144]=-0.0905455864075847;
  coords[145]=0.3000939329402901;
  coords[146]=-0.22837926412108878;
  coords[147]=-0.20711334723323238;
  coords[148]=0.45017480408390376;
  coords[149]=-0.49492932461998207;
  coords[150]=0.041428457629927085;
  coords[151]=-0.48833001850612767;
  coords[152]=0.1018317078007771;
  coords[153]=-0.48768343494393396;
  coords[154]=0.36961018266280654;
  coords[155]=0.4773305012741377;
  coords[156]=0.22255818378285508;
  coords[157]=0.36750039563306314;
  coords[158]=-0.48421438420511276;
  coords[159]=0.14953479103762057;
  coords[160]=-0.1484083780302572;
  coords[161]=-0.3520389116514898;
  coords[162]=-0.04918611662845995;
  coords[163]=-0.4515623105888268;
  coords[164]=0.008675167218021973;
  coords[165]=0.30502459042699503;
  coords[166]=0.0691383128891293;
  coords[167]=0.09106225239469867;
  coords[168]=0.08403241809528572;
  coords[169]=-0.11580566981396712;
  coords[170]=0.4337974740651488;
  coords[171]=0.4831195274103003;
  coords[172]=-0.02419354659538453;
  coords[173]=-0.2626185851094771;
  coords[174]=-0.19281536886551864;
  coords[175]=-0.3694676725734728;
  coords[176]=-0.26170104262540295;
  coords[177]=-0.4819015900272876;
  coords[178]=0.28020829951686366;
  coords[179]=-0.21019039676593632;
  coords[180]=-0.23295845777744606;
  coords[181]=0.33354885063650286;
  coords[182]=0.02832680855630998;
  coords[183]=0.3375508104939249;
  coords[184]=-0.03432461464308212;
  coords[185]=0.416802163696285;
  coords[186]=0.09499991399934982;
  coords[187]=-0.39325708950085647;
  coords[188]=-0.08518805988705602;
  coords[189]=0.0777512383251161;
  coords[190]=0.13142135411586453;
  coords[191]=0.21114495877922645;
  coords[192]=-0.2055082929781975;
  coords[193]=-0.006691582021893239;
  coords[194]=0.20596565610682815;
  coords[195]=0.15893836716919474;
  coords[196]=-0.30044613929233965;
  coords[197]=-0.39469751989696567;
  coords[198]=-0.3623960831443602;
  coords[199]=-0.15704651304619832;
  qsort(coords, 100, 2*sizeof(double), cmp);

  for (int i = 0; i < 100; ++i)
    printf("%8g %8g\n", coords[2*i], coords[2*i+1]);
  free(coords);
  return 0;
}
#endif
