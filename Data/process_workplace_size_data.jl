#=
Purpose:
Process data from the ONS 'UK business: activity, size and location' dataset
Produce breakdowns of the fraction of workplaces within specified employee size ranges
Separate estimates produced for different work sectors
=#
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# SET PATHS & LOAD ENVIRONMENT
#-------------------------------------------------------------------------------

#Set paths
cd(dirname(@__FILE__))

#Load environment
using Pkg
Pkg.activate("../")

#-------------------------------------------------------------------------------
# LOAD PACKAGES
#-------------------------------------------------------------------------------
using XLSX
using DelimitedFiles

#-------------------------------------------------------------------------------
# LOAD DATA
#-------------------------------------------------------------------------------
work_division_data = XLSX.readdata("ukbusinessworkbook2020.xlsx", "Table 3", "B8:H95")
work_classes_data = XLSX.readdata("ukbusinessworkbook2020.xlsx", "Table 4", "B7:H621")

#-------------------------------------------------------------------------------
# SET UP MAPPINGS TO THE 41 WORK SECTORS
#-------------------------------------------------------------------------------
n_bins = size(work_division_data,2)
work_sector_data = zeros(Int64,41,n_bins)

# 1. Agriculture                12. Postal                              23. Employment and HR                   34. Betting
# 2. Mining                     13. Accommodation                       24. Travel Agency                       35. Sport
# 3. Manufacturing (food)       14. Restaurant and Bar                  25. Security                            36. Theme Parks
# 4. Manufacturing (other)      15. Broadcasting and Communications     26. Cleaning                            37. Religious and Political Organisations
# 5. Utilities and Waste        16. Information Technology              27. Office (other)                      38. Repair
# 6. Construction               17. News                                28. Public Administration and Defence   39. Hairdressers
# 7. Motor Trade                18. Banking/Accounting                  29. Education                           40. Funeral
# 8. Wholesale                  19. Real Estate                         30. Hospital/Doctor/Dental              41. Personal Services
# 9. Retail                     20. Professional/Science/Tech           31. Care Homes
# 10. Transport                 21. Veterinary                          32. Social Work
# 11. Transport Support         22. Rental Companies                    33. Arts
work_sectors_from_division_data_idxs = [collect(1:34);37;38]
start_idx = [1,4,9,11,32,38,41,42,43,44,47,48,49,50,51,54,56,57,60,61,67,68,69,70,71,72,73,74,75,76,77,78,79,81,83,84]
end_idx = [3,8,10,31,37,40,41,42,43,46,47,48,49,50,53,55,56,59,60,66,67,68,69,70,71,72,73,74,75,76,77,78,80,81,83,84]
for work_sector_itr = 1:length(work_sectors_from_division_data_idxs)
    row_to_assign_to_idx = work_sectors_from_division_data_idxs[work_sector_itr]
    work_sector_data[row_to_assign_to_idx,:] = sum(work_division_data[start_idx[work_sector_itr]:end_idx[work_sector_itr],:], dims=1)
end

# Assign sports, theme parks, hairdressers, funeral data, personal services
# from work_classes_data
work_sector_data[35,:] = sum(work_classes_data[587:590,:], dims=1)
work_sector_data[36,:] = sum(work_classes_data[591:592,:], dims=1)
work_sector_data[39,:] = sum(work_classes_data[608:608,:], dims=1)
work_sector_data[40,:] = sum(work_classes_data[609:609,:], dims=1)
work_sector_data[41,:] = sum(work_classes_data[[607;collect(610:615)],:], dims=1)

#-------------------------------------------------------------------------------
# PRODUCE CDF VALUES
#-------------------------------------------------------------------------------

work_sector_data_cdfs = cumsum(work_sector_data, dims=2)./sum(work_sector_data, dims=2)

#-------------------------------------------------------------------------------
# SAVE ARRAY TO FILE
#-------------------------------------------------------------------------------
#Take tranpose to give row by workplace size and column by work sector
writedlm("../../Data/worker_model/work_sector_proportions_by_size.csv",work_sector_data_cdfs')
