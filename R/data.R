#' Patient Admission data
#'
#' @source Data from a set of hospitals in Victoria (Australia)
#' @format A data frame with columns:
#' \describe{
#'  \item{LOS.total}{Lenght of stay}
#'  \item{Age}{Age of patient}
#'  \item{NumberEpisodes}{Number of hospital episodes}
#'  \item{RealAdmissionDate}{Date of admission}
#'  \item{RealSeparationDate}{Date of discharge}
#'  \item{GenderName}{Gender of patient}
#'  \item{Postcode}{Postcode of patient}
#'  \item{MaritalStatusNameCorr}{Marital status of patient}
#'  \item{DiseaseTypesCorr}{Disease type - 9 levels}
#' }
#' @examples
#' \dontrun{
#'  PatientAdmission
#' }
"PatientAdmission"

#' Patient medication data
#' @format A data frame with columns:
#' \describe{
#'  \item{ofp}{}
#'  \item{ofnp}{}
#'  \item{opp}{}
#'  \item{opnp}{}
#'  \item{emer}{}
#'  \item{hosp}{}
#'  \item{health}{}
#'  \item{numchron}{}
#'  \item{adldiff}{}
#'  \item{region}{}
#'  \item{age}{}
#'  \item{black}{}
#'  \item{gender}{}
#'  \item{married}{}
#'  \item{school}{}
#'  \item{faminc}{}
#'  \item{employed}{}
#'  \item{privins}{}
#'  \item{medicaid}{}
#' }
#' @examples
#' \dontrun{
#'  DebTrivedi
#' }
"DebTrivedi"