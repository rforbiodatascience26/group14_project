

# The following script will process the metadata files for analysis in:
#    - Global proteomics
#    - Phosphoproteomics
#    - Acetyl proteomics

# There is a lot of extra variables whose effect we don't want to test on this study

# The variable of study will be:
#     a) Endurance training, samples taken at the end of an 8 week endurance training

# Against the control (untrained)

# Confounding variables to take into account:
#     b) Sex
#     c) Batch (TmT plex run)


# Section 1:  [ Select vial IDs for study ]:

# We only want to include in our study the samples of mouses that have undergone
# an endudrance training of 8 weeks and of mouses who have not trained during that period of time

# To do that, we must select the following values for the variables :

# To select control mice:
    # key___intervention == 2
    # study_group_timepoint == "Control - 8 weeks"

# To select trained mice:
    # key___intervention == 3
    # study_group_timepoint == "Training - 8 weeks"

clean_metadata_table <- meta_vial_label |> 
  filter((
    key___intervention == 3 & study_group_timepoint == "Control - 8 weeks")
    |  key___intervention == 2 & study_group_timepoint == "Training - 8 weeks") |> 
  
  # Section 2 : [ Generate the new metadata table ]:
  
  select(viallabel,
         key___intervention,    # 2 = Endurance training 3 = control
         registration___sex)|>  # 1= male 2= female
  
  mutate(
    condition = case_when(
    key___intervention == 2 ~ "Control",
    key___intervention == 3 ~ "Training"),
    
    sex = case_when(
    registration___sex == 1 ~ "F", 
    registration___sex == 2 ~ "M")) |> 
  
  select(viallabel,condition, sex)
  
