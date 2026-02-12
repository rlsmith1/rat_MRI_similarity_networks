
#==============================================================================#
# Find the center of mass for each Waxholm ROI
# This script uses AFNI command line functions
#==============================================================================#


# Setup -------------------------------------------------------------------


library(tidyverse)
base_dir <- "~/Documents/PhD/projects/CamRat/CamRat/"
objects_dir <- paste0(base_dir, "outputs/objects")

atlas_dir <- "~/Documents/PhD/projects/CamRat/atlases/WHS_SD_rat_atlas_v4_pack/"
df_labels <- read.table(paste0(atlas_dir, "WHS_SD_rat_atlas_v4.label")) %>% 
  dplyr::rename_all(~c("IDX", "-R-",  "-G-", "-B-", "-A--",  "VIS", "MSH", "region_of_interest")) %>% 
  as_tibble()


## Get IDX values from WHS
df_labels %>% 
  pull(IDX) %>% 
  unique


# AFNI command line -------------------------------------------------------



## Separate right and left hemispheres
# 3dcalc -a WHS_SD_rat_atlas_v4.nii -RAI -expr 'a*step(x)' -prefix left_mask.nii (left hemisphere only)
# 3dcalc -a WHS_SD_rat_atlas_v4.nii -RAI -expr 'a*step(-x)' -prefix right_mask.nii (right hemisphere only)

# 3dCM -Icent -local_ijk -roi_vals 1 3 4 5 6 7 10 32 33 34 35 36 37 38 40 41 42 43 45 46 47 48 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 93 94 95 96 97 98 99 100 108 109 110 112 113 114 115 119 120 121 122 123 125 126 127 128 129 130 131 132 133 134 135 136 137 138 139 140 141 142 143 145 146 150 151 152 153 157 158 159 160 162 163 164 180 181 182 183 184 187 188 189 192 193 195 196 197 198 199 200 201 204 205 206 207 208 210 211 213 214 215 216 218 219 221 222 223 227 228 229 230 231 232 233 235 236 238 239 240 242 246 247 248 249 254 255 257 260 266 267 268 270 272 278 280 281 282 283 284 285 286 287 290 291 293 294 295 297 298 299 400 401 402 403 404 405 406 407 408 409 410 411 412 413 414 416 417 418 420 422 423 424 425 427 429 430 432 433 436 442 443 444 448 500 501 502 left_mask.nii > 3dCM_local_ijk_LH.out



# Read results ------------------------------------------------------------


# Read results
lines_LH <- readLines(paste0(base_dir, "data/WHS_analysis/3dCM_local_ijk_LH.txt"))
lines_LH <- lines_LH[!grepl("^#Dset", lines_LH)]

lines_RH <- readLines(paste0(base_dir, "data/WHS_analysis/3dCM_local_ijk_RH.txt"))
lines_RH <- lines_RH[!grepl("^#Dset", lines_RH)]


# Combine left and right hemispheres into tibble
df_roi_centers <- tibble(
  hemisphere = "left",
  IDX = lines_LH[str_detect(lines_LH, "^#ROI")] %>% str_remove("#ROI ") %>% as.integer,
  center = lines_LH[!str_detect(lines_LH, "^#ROI")]
) %>% 
  bind_rows(
    tibble(
      hemisphere = "right",
      IDX = lines_RH[str_detect(lines_RH, "^#ROI")] %>% str_remove("#ROI ") %>% as.integer,
      center = lines_RH[!str_detect(lines_RH, "^#ROI")]
    )
  ) %>% 
  separate(center, into = c("x", "y", "z"), sep = "\\s+") %>% 
  mutate_at(vars(c(x, y, z)), ~ as.integer(.x)) %>% 
  left_join(df_labels %>% 
              dplyr::select(IDX, region_of_interest)
  ) %>% 
  mutate(region_of_interest = str_trim(region_of_interest))

# add ROI label information
save(df_roi_centers, file = paste0(objects_dir, "WHS_centers.RDS"))


