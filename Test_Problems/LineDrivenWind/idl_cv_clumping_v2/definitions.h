#define  PHYSICS                        HD
#define  DIMENSIONS                     2
#define  GEOMETRY                       SPHERICAL
#define  BODY_FORCE                     VECTOR
#define  COOLING                        BLONDIN
#define  RECONSTRUCTION                 PARABOLIC
#define  TIME_STEPPING                  RK3
#define  NTRACER                        1
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            15
#define  LINE_DRIVEN_WIND               SIROCCO_MODE

/* -- physics dependent declarations -- */

#define  DUST_FLUID                     NO
#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 ALWAYS
#define  INCLUDE_LES                    NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  RADIATION                      NO
#define  ROTATING_FRAME                 NO

/* -- user-defined parameters (labels) -- */

#define  MU                             0
#define  RHO_0                          1
#define  R_0                            2
#define  RHO_ALPHA                      3
#define  CENT_MASS                      4
#define  DISK_MDOT                      5
#define  T_ISO                          6
#define  L_star                         7
#define  f_x                            8
#define  f_uv                           9
#define  T_x                            10
#define  KRAD                           11
#define  ALPHARAD                       12
#define  DFLOOR                         13
#define  GAMMA                          14

/* [Beg] user-defined constants (do not change this line) */

#define  WARNING_MESSAGES               NO
#define  INTERNAL_BOUNDARY              YES
#define  INTERNAL_BOUNDARY_CFL          YES
#define  SHOCK_FLATTENING               MULTID
#define  CHAR_LIMITING                  YES
#define  LIMITER                        VANLEER_LIM
#define  FAILSAFE                       YES
#define  PPM_ORDER                      4
#define  UNIT_DENSITY                   1e-9
#define  UNIT_LENGTH                    8.7e8
#define  UNIT_VELOCITY                  3.1e8
#define  UNIT_MASS                      (UNIT_DENSITY*UNIT_LENGTH*UNIT_LENGTH*UNIT_LENGTH)
#define  UNIT_ACCELERATION              (UNIT_VELOCITY*UNIT_VELOCITY/UNIT_LENGTH)
#define  UNIT_FORCE                     (UNIT_MASS*UNIT_ACCELERATION)
#define  UNIT_TIME                      (UNIT_LENGTH/UNIT_VELOCITY)
#define  UNIT_PRESSURE                  (UNIT_DENSITY*UNIT_VELOCITY*UNIT_VELOCITY)

/* [End] user-defined constants (do not change this line) */
