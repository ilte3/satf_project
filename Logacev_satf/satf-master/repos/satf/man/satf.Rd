\name{em2}
\alias{em2}
\title{ Fit SATF function to response-signal data. }
\description{
  Fit SATF function to response-signal data. 
}


\usage{
  satf(dv, signal, start, contrasts, bias, data, time, metric,                   
     constraints = list(), trial.id = NULL, reoptimize.times = 1, method = "SUBPLEX", debug = F,
     optimize.incrementally = FALSE, reoptimize.dprime = TRUE, 
     reoptimize.criterion = TRUE, reoptimize.corr = TRUE)

(m.raw <- satf(dv=c(response=~response), signal=~signal, 
               start = c(asymptote=2, invrate=1, intercept=.4, bias.max=1, bias.invrate=1, bias.intercept=-1, bias.min=0),
               contrasts = c(asymptote=~1+c2, invrate=~1+c2, intercept=~1+c2),
               bias = ~1+c2, constraints=list(), time=~time, trial.id=~trial.id,
               data=data, metric="logLikRaw", debug=T, method="Nelder-Mead", optimize.incrementally=FALSE, reoptimize.times=5))

}
\arguments{
  \item{dv}{ A named list with one or two formulae specifying the dependent variable column(s). Either list(response=~...), or list(n.responses.yes=~..., n.responses=~...).
  %, or list(dprime=...).
  }
  \item{signal}{A formula specifying a column which indicates whether the datapoint is from a signal trial or from a noise trial.}
  \item{contrasts}{A named list specifying the contrasts for SATF parameters. Contrasts for 'asymptote', 'invrate', and 'intercept' have to be provided. For example, contrasts=list(asymptote=~1, invrate=~1, intercept=~1) requests only one common asymptote to be estimated.}
  \item{start}{A named vector with start values, e.g. c(asymptote=2, invrate=1, intercept=.5). }

  \item{constraints}{A named list of constraints for all coeficients. Default constraints set the asymptote, invrate and intercept parameters to be strictly positive. Constraints can apply to one parameters (e.g., constraints=list('asymptote'=c(0, 2))), or can use regular expressions to constrain several parameters (e.g., constraints=list('asymptote.c.*'=c(0, 2))).
  }
  \item{time}{A formula specifying the time column (e.g., time=~time).}
  \item{trial.id}{A formula specifying the trial id column  (e.g., trial.id=~trial). The trial id is necessary to account for within-trial correlation of responses if the data is from a MR-SAT experiment.}

  \item{bias}{Either a formula, or a named list specifying contrasts for the decition criterion. If a named list is provided, the components (bias.max, bias.invrate, bias.intercept, bias.min) have to provide formulae specifying the constrasts for these parameters. If the contrasts are the same for all parameters of the criterion function, one formula can be provided. (The bias is ) }
  \item{data}{data.frame containing the data.}
  \item{metric}{String specifying the metric by which to fit the model. }

  \item{plot}{NULL, or list providing the details necessary for plotting the data.}
  \item{metric}{}

  \item{fixationdur}{A vector of fixation durations, typically in milliseconds. The elements of the vector are considered integers.}

  \item{trialinfo}{A data.frame containing trial-level information that shall be incorporated into the resulting data.frame with the eyetracking measures. Typically, it will consist of subject id, item id, condition and group/list. Information about regions cannot be included at this point, but possibly in future releases.}

  \item{measures}{The type of measures you want to compute. Possible values are "standard" or "CRI". "standard" will return the common measures used in psycholinguistics, and some not-so-common ones. "CRI" will return "conditional regression information", see the 'Details' section for details.}

  \item{regressiveFirstPass}{If 'regressiveFirstPass' is set to FALSE (default behavior), the first-pass reading time and associated measures will be 0 for regions which were fixated for the first time *after* a fixation on a region to the right (i.e., a region with a larger 'rid').
  If 'regressiveFirstPass' is set to TRUE, such first fixations are treated like ordinary first fixations.}
		    
  	\item{assimilateSaccades}{If TRUE, em() adds saccade durations to fixation durations. Fixation start and end times have to be provided.}
    
	\item{fixationstart, fixationend}{Vectors of integers containing start and end times of  fixations. They are not used unless 'assimilateSaccades' is set to TRUE.}
  }
\details{

The em() function computes psycholinguistic measures used in research on sentence processing. 

When the 'measures' parameter is set to "standard" the following measures are computed:
first fixation duration (FFD), single fixation duration (SFD), first pass reading time or gaze duration (FPRT), right-bounded reading time (RBRT), regression path duration or go-past duration (RPD), total fixation time (TFT), rereading time (RRT), rereading time regressive (RRTR), rereading time progressive (RRTP), first fixation progressive (FFP), total regression count (TRC), right-bounded regression count (RBRC). 
When 'measures' is set to "CRI" then "conditional regression information" will be returned. CRI consists of two components: CFC (conditional fixation count) and CFT (conditional fixation time). CRI provides detailed information about the targets of regressive fixations which started at certain positions. The CFT for a 'source region' and a 'destination region' provides the total time spent fixating the destination region *after* a fixation on the source, but *before* any fixation on a position to the right of the source. CFC provides the number of fixations on the destination region during this perion of time.

}
\value{
When measures="standard", the following measures will be computed for each trial and each region of interest (rid):
 
   \item{FFD}{(First Fixation Duration)
	 Duration of the first fixation on a position if (and only if) the fixation was progressive.
	 Zero otherwise.}

   \item{SFD}{(Single Fixation Duration)
	Duration of the fixation on a position if it was the *only* fixation on this region, i.e. if no subsequent fixation on this position followed. Zero if there were several fixations on this region.}

   \item{FPRT}{(First Pass Reading Time, Gaze Duration)
   	Sum of all first-pass fixation durations on a region before *any* other region is fixated. (What exactly constitutes a first pass is determined by the parameter 'regressiveFirstPass'.)}

   \item{RBRT}{(Right Bounded Reading Time)
   	Sum of all first-pass fixation durations on a position before another position to the *right* is fixated. (What exactly constitutes a first pass is determined by the parameter 'regressiveFirstPass'.)}

   \item{RPD}{(Regression Path Duration, Go-Past Duration)
   	Sum of all first-pass fixation durations on a position n and all preceding positions in the time period between the first fixation on n and the first fixation on anything to the right of n.}
   
   \item{CRPD}{(Cumulative Regression Path Duration)
	The CRPD of position n is the total amount of time a participant spent reading the sentence until reaching region n+1. It is the sum of the RPDs of all preceding regions and the RPD of the current region. }   

   \item{TFT}{(Total Fixation Time)
   	Sum of all fixation durations on a region.}

   \item{RRT}{(Re-reading Time)
	Sum of all second-pass fixation durations. (RRT = TFT - FPRT)}

   \item{RRTR}{(ReReading Time Regressive)
	Sum of all second-pass fixation durations on a position that occured *after* a fixation on a region further to the right. (RRTR = TFT - RBRT)}

   \item{RRTP}{(ReReading Time Progressive)
	Sum of all second-pass fixation durations on a position that took place *before* a fixation on a region further to the right. (RRTP = RBRT - FPRT) }

   \item{FFP}{(First Fixation Progressive)
	0 if material downstream was viewed before the first fixation on this position, 
	1 otherwise. }

   \item{TRC}{ (Total Regression Count)
	Total number of regressions from this position.}

   \item{RBRC}{(Right-Bounded Regression Count)
	Number of regressions from this position given *before* any region further to the right has been fixated.}
 
When measures="CRI", the following measures are computed for each trial and each pair of destination region (rid) and source source region (sourcerid):

   \item{CFC}{(Conditional Fixation Count)
	Number of fixations on rid after sourcerid has been fixated, but before any region right of sourcerid has been fixated. Two subsequent fixations on one region (i.e. when no other region is looked at inbetween) are considered one fixation.}
   \item{CFT}{(Conditional Fixation Time)
	Amount of time spent on fixating rid under the conditions described for CFC. }
 
 
}
\references{ 
Rayner, K. (1998). Eye movements in reading and information processing: 20 years of research. Psychological bulletin, 124(3), 372.

Clifton, C., Staub, A., & Rayner, K. (2007). Eye movements in reading words and sentences. Eye movements: A window on mind and brain, 341-372.

Mitchell, D. C., Shen, X., Green, M. J., & Hodgson, T. L. (2008). Accounting for regressive eye-movements in models of sentence processing: A reappraisal of the Selective Reanalysis hypothesis. Journal of Memory and Language, 59(3), 266-293.
}
\author{ Pavel Logacev <pavel.logacev@gmail.com>, Shravan Vasishth <vasishth@acm.org> }
%\seealso{ \code{\link{summary.maxchi}}, \code{\link{phylpro}}  }
\examples{ 
# generate a toy dataset
data <- em2.dset()

# generate standard measures
data.em <- em2(data$rid, data$fixationdur, data[,c('subject','trial','response')])

# generate CRI 
data.cri <- em2(data$rid, data$fixationdur, data[,c('subject','trial','response')], measures='CRI')
	
}
%\keyword{ htest }
