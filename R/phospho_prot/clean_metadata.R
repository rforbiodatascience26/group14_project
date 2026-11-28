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
      registration___sex == 2 ~ "M"),
    
    vial_label = as.character(viallabel)) |> 
  
  select(vial_label,condition, sex)


table(clean_metadata_table$condition)

table(clean_metadata_table$sex)

filter(clean_metadata_table, vial_label == "90571015305")
