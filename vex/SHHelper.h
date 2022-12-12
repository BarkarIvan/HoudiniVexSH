#ifndef ME_SPHERICAL_HARMONICS
#define ME_SPHERICAL_HARMONICS

int factorial(int i) {
    int factorial_table[] = array(
    1,
    1,
    2,
    6,
    24,
    120,
    720,
    5040,
    40320,
    362880,
    3628800,
    39916800,
    479001600,
    6227020800,
    87178291200,
    1307674368000,
    20922789888000,
    355687428096000,
    6402373705728000,
    121645100408832000,
    2432902008176640000,
    51090942171709440000,
    1124000727777607680000,
    25852016738884976640000,
    620448401733239439360000,
    15511210043330985984000000,
    403291461126605635584000000,
    10888869450418352160768000000,
    304888344611713860501504000000,
    8841761993739701954543616000000,
    265252859812191058636308480000000,
    8222838654177922817725562880000000,
    263130836933693530167218012160000000,
    8683317618811886495518194401280000000
    );

    return factorial_table[i];
}

float P(int l; int m; float x)
{
  // evaluate an Associated Legendre Polynomial P(l,m,x) at x
  float pmm = 1.0;
  if(m>0) {
    float somx2 = sqrt((1.0-x)*(1.0+x));
    float fact = 1.0;
    int i;
    for(i=1; i<=m; i++) {
      pmm *= (-fact) * somx2;
      fact += 2.0;
    }
  }
  if(l==m) return pmm;
  float pmmp1 = x * (2.0*m+1.0) * pmm;
  if(l==m+1) return pmmp1;
  float pll = 0.0;
  int ll;
  for(ll=m+2; ll<=l; ++ll) {
    pll = ( ((float)(2.0*ll) - 1.0)*x*pmmp1 - ((float)(ll+m) - 1.0)*pmm ) / (ll - m);
    pmm = pmmp1;
    pmmp1 = pll;
  }
  return pll;
}


float K(int l; int m)
{
  // renormalisation constant for SH function
  float temp = ((2.0*l+1.0)*factorial(l-m)) / (4.0*PI*factorial(l+m));
  return sqrt(temp);
}

float SH(int l; int m; float theta; float phi)
{
  // return a point sample of a Spherical Harmonic basis function
  // l is the band, range [0..N]
  // m in the range [-l..l]
  // theta in the range [0..Pi]
  // phi in the range [0..2*Pi]

  float sqrt2 = sqrt(2.0);

  if(m==0)
    return K(l,0)*P(l,m,cos(theta));
  else if(m>0)
    return sqrt2*K(l,m)*cos(m*phi)*P(l,m,cos(theta));
  else
    return sqrt2*K(l,-m)*sin(-m*phi)*P(l,-m,cos(theta));
}

int getindex(int l; int m) {
    return (l*(l+1))+m;
}

float mtx_get(matrix3 matx; int index) {

    int j, k;
    // hard coded for 3x3 matrix
    j = index / 3;
    k = index % 3;

    return getcomp(matx, j, k);
}

void mtx_set(matrix3 mtx; int index; float val) {
    int j, k;

    // hard coded for 3x3 matrix
    j = index / 3;
    k = index % 3;

    setcomp(mtx, val, j, k);
}


//spherical {r, θ, φ} to cartesian {x, y, z}
function vector toCartesian(vector spherical;)
{
    float r = spherical.x;
    float theta = spherical.y;
    float phi = spherical.z;

    float x = r * sin(theta) * cos(phi);
    float y = r * cos(theta);
    float z = r * sin(theta) * sin(phi);

    return set(x, y, z);
}


//cartesian {x, y, z} to spherical {r, θ, φ}
function vector toSpherical(vector cartesian;)
{
    float r = length(cartesian); // r = √(x² + y² + z²)
    float theta = acos(cartesian.y/r);
    float phi = atan2(cartesian.z,cartesian.x); //use atan2 for built-in checks


    return set(r, theta, phi);
}


function vector UVtoPolar(vector pos)
{
  
  float radius = sqrt(pos.x*pos.x + pos.y*pos.y + pos.z*pos.z);
  float u = atan2(pos.y, pos.x) + PI;
  float v = acos(pos.z / radius);
  u = fit(u, 0, 2*PI, 0, 1); 
  v = v/ PI;
  return set(u, v, radius);	
}

#endif /* ME_SPHERICAL_HARMONICS */
