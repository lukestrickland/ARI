
# Constants for experiment types
TASKTYPE_ARI            <- 1
TASKTYPE_FLASHING       <- 2
TASKTYPE_FLASHING_ELISE <- 3

# Constants for trial types
# AR/ARI
TASKCODE_ARI_AR_TRIAL  <- 1
TASKCODE_ARI_ARI_TRIAL <- 2

# FLASHING DOTS
TASKCODE_FLASHING_L_EASY_TRIAL <- 1
TASKCODE_FLASHING_R_EASY_TRIAL <- 2
TASKCODE_FLASHING_L_HARD_TRIAL <- 3
TASKCODE_FLASHING_R_HARD_TRIAL <- 4
TASKCODE_FLASHING_STOP_TRIAL   <- 8   # NOTE: This is bitwise ORed with previous settings

# Constants for block types
BLOCKTYPE_AR_ONLY       <- 1
BLOCKTYPE_AR_ARI        <- 2
BLOCKTYPE_FLASHING_ONLY <- 3
BLOCKTYPE_FLASHING_SS   <- 4

# For Elise's experiment
BLOCKTYPE_FLASHING_CONSTANT_POINTS   <- 10
BLOCKTYPE_FLASHING_CONSTANT_MONEY    <- 11
BLOCKTYPE_FLASHING_JACKPOT_POINTS    <- 12
BLOCKTYPE_FLASHING_JACKPOT_MONEY     <- 13
BLOCKTYPE_FLASHING_PROPORTIONAL      <- 14

BLOCKTYPE_SHOW_INSTRUCTIONS <- 64           # Luke's experiment / NOTE: This is bitwise ORed with previous settings

BLOCKTYPE_SHOW_INSTRUCTIONS_PRACTICE <- 128 # Elise's experiment / NOTE: This is bitwise ORed with previous settings
BLOCKTYPE_SHOW_INSTRUCTIONS_TEST1    <- 256 # Elise's experiment / NOTE: This is bitwise ORed with previous settings
BLOCKTYPE_SHOW_INSTRUCTIONS_TEST2    <- 512 # Elise's experiment / NOTE: This is bitwise ORed with previous settings
