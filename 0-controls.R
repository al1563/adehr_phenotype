
library(DBI)
library(dplyr)
library(odbc)
library(MatchIt)
library(ggplot2)

con <- dbConnect(odbc::odbc(), .connection_string = 'REDACTED')

threshold_encounters <- FALSE

# Query AD cohort
alz_cohort <- dbGetQuery(con, 'INSERT QUERY HERE')
alz_cohort[,'isad'] <- rep('1',nrow(alz_cohort))
sprintf("alz_cohort dimensions: %s %s", dim(alz_cohort)[1],dim(alz_cohort)[2])
head(alz_cohort)

# Query all other patients
background_cohort <- dbGetQuery(con, 'INSERT QUERY HERE')
background_cohort[,'isad'] <- rep('0',nrow(alz_cohort))
sprintf("background_cohort dimensions: %s %s", dim(background_cohort)[1],dim(background_cohort)[2])

# put them in one giant dataframe. isaad column labels 1 as alzheimer cohort, and 0 as background cohort.
alz_background_pts <- rbind(alz_cohort,background_cohort)
alz_background_pts <- alz_background_pts[complete.cases(alz_background_pts),]
print(sum(is.na(alz_background_pts))) # Make sure no missing values

# run MatchIt for propensity score matching.
start_time <- Sys.time()

if(threshold_encounters){
    m.out = matchit(isalz ~ FirstRace + estimated_age + Sex + Status + encCnt + dateDif, 
      data = alz_background_pts, method = "nearest", ratio = 2)
    m.data <- match.data(m.out,distance='pscore')
    
    print('Saving plot...')
    pdf(paste0('Data/encPSplot.pdf'))
    myplot <- ggplot(m.data, aes(x=pscore, fill=as.factor(isalz))) + geom_density(alpha=0.25)+ xlab("Estimated PS")
    print(myplot)

    print('Saving data...')
    save(m.data, m.out, file = paste0('Data/controlencdata.RData'))
    write.csv(m.data,paste0('Data/controlencdata.csv'))
} else{
    m.out = matchit(isalz ~ FirstRace + estimated_age + Sex + Status, 
      data = alz_background_pts, method = "nearest", ratio = 2)
    m.data <- match.data(m.out,distance='pscore')
    print('Saving plot...')
    pdf(paste0('Data/PSplot.pdf'))
    myplot <- ggplot(m.data, aes(x=pscore, fill=as.factor(isalz))) + geom_density(alpha=0.25)+ xlab("Estimated PS")
    print(myplot)

    print('Saving data...')
    save(m.data, m.out, file = paste0('Data/controldata.RData'))
    write.csv(m.data,paste0('Data/controldata.csv'))
}

end_time <- Sys.time()
print(end_time - start_time)

