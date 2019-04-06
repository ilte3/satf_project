\name{Auditory_demo}
\alias{Auditory_demo}
\docType{data}
\title{
A data from single subject undertaking auditory version of MR-SAT experiment
}
\description{
A data ...
}
\usage{data(Auditory_demo)}
\format{
  A data frame with 480 observations on the following 99 variables.
  \describe{
    \item{\code{ExperimentName}}{a character vector. The name of the experiment as specified in the Eprime script.}
    \item{\code{Subject}}{a numeric vector}
    \item{\code{Session}}{a numeric vector}
    \item{\code{Display.RefreshRate}}{a numeric vector. The refresh rate of the monitor used during the experiment}
    \item{\code{Group}}{a numeric vector}
    \item{\code{RandomSeed}}{a numeric vector. The random seeds used during the experiment.}
    \item{\code{SubjectID}}{a numeric vector. The ID of the participant, as specified at the beginning of the experiement.}
    \item{\code{Block}}{a numeric vector}
    \item{\code{Condition}}{a numeric vector}
    \item{\code{CorrResp}}{a numeric vector. Correct key }
    \item{\code{Item}}{a numeric vector}
    \item{\code{key01}}{The 1st key pressed.}
    \item{\code{key02}}{The 2nd key pressed.}
    \item{\code{key03}}{The 3rd key pressed.}
    \item{\code{key04}}{The 4th key pressed.}
    \item{\code{key05}}{The 5th key pressed.}
    \item{\code{key06}}{The 6th key pressed.}
    \item{\code{key07}}{The 7th key pressed.}
    \item{\code{key08}}{The 8th key pressed.}
    \item{\code{key09}}{The 9th key pressed.}
    \item{\code{key10}}{The 10th key pressed.}
    \item{\code{key11}}{The 11th key pressed.}
    \item{\code{key12}}{The 12th key pressed.}
    \item{\code{key13}}{The 13th key pressed.}
    \item{\code{key14}}{The 14th key pressed.}
    \item{\code{key15}}{The 15th key pressed.}
    \item{\code{key16}}{The 16th key pressed.}
    \item{\code{key17}}{The 17th key pressed.}
    \item{\code{key18}}{The 18th key pressed.}
    \item{\code{key19}}{The 19th key pressed.}
    \item{\code{key20}}{The 20th key pressed.}
    \item{\code{key21}}{The 21st key pressed.}
    \item{\code{key22}}{The 22nd key pressed.}
    \item{\code{key23}}{The 23rd key pressed.}
    \item{\code{key24}}{The 24th key pressed.}
    \item{\code{key25}}{The 25th key pressed.}
    \item{\code{key26}}{The 26th key pressed.}
    \item{\code{key27}}{The 27th key pressed.}
    \item{\code{key28}}{The 28th key pressed.}
    \item{\code{key29}}{The 29th key pressed.}
    \item{\code{key30}}{The 30th key pressed.}
    \item{\code{LastWord}}{a character vector}
    \item{\code{List1}}{a numeric vector}
    \item{\code{List1.Cycle}}{a numeric vector}
    \item{\code{List1.Sample}}{a numeric vector}
    \item{\code{NotBoth}}{a character vector}
    \item{\code{OrientSound.OffsetTime}}{a numeric vector}
    \item{\code{OrientSound.OnsetTime}}{a numeric vector}
    \item{\code{Procedure}}{a character vector}
    \item{\code{ResponseCollector.OffsetTime}}{a numeric vector}
    \item{\code{ResponseCollector.OnsetTime}}{a numeric vector}
    \item{\code{ResponseCollector.StartTime}}{a numeric vector}
    \item{\code{rt01}}{The response latency of the 1st key press.}
    \item{\code{rt02}}{The response latency of the 2nd key press.}
    \item{\code{rt03}}{The response latency of the 3rd key press.}
    \item{\code{rt04}}{The response latency of the 4th key press.}
    \item{\code{rt05}}{The response latency of the 5th key press.}
    \item{\code{rt06}}{The response latency of the 6th key press.}
    \item{\code{rt07}}{The response latency of the 7th key press.}
    \item{\code{rt08}}{The response latency of the 8th key press.}
    \item{\code{rt09}}{The response latency of the 9th key press.}
    \item{\code{rt10}}{The response latency of the 10th key press.}
    \item{\code{rt11}}{The response latency of the 11th key press.}
    \item{\code{rt12}}{The response latency of the 12th key press.}
    \item{\code{rt13}}{The response latency of the 13th key press.}
    \item{\code{rt14}}{The response latency of the 14th key press.}
    \item{\code{rt15}}{The response latency of the 15th key press.}
    \item{\code{rt16}}{The response latency of the 16th key press.}
    \item{\code{rt17}}{The response latency of the 17th key press.}
    \item{\code{rt18}}{The response latency of the 18th key press.}
    \item{\code{rt19}}{The response latency of the 19th key press.}
    \item{\code{rt20}}{The response latency of the 20th key press.}
    \item{\code{rt21}}{The response latency of the 21st key press.}
    \item{\code{rt22}}{The response latency of the 22nd key press.}
    \item{\code{rt23}}{The response latency of the 23rd key press.}
    \item{\code{rt24}}{The response latency of the 24th key press.}
    \item{\code{rt25}}{The response latency of the 25th key press.}
    \item{\code{rt26}}{The response latency of the 26th key press.}
    \item{\code{rt27}}{The response latency of the 27th key press.}
    \item{\code{rt28}}{The response latency of the 28th key press.}
    \item{\code{rt29}}{The response latency of the 29th key press.}
    \item{\code{rt30}}{The response latency of the 30th key press.}
    \item{\code{Running}}{a character vector}
    \item{\code{Sentence}}{The relative path to the audio recording of the stimulus sentences.}
    \item{\code{SentencePlayer.ActionDelay}}{a numeric vector}
    \item{\code{SentencePlayer.Duration}}{a numeric vector}
    \item{\code{SentencePlayer.FinishTime}}{a numeric vector}
    \item{\code{SentencePlayer.OffsetTime}}{a numeric vector}
    \item{\code{SentencePlayer.OnsetTime}}{a numeric vector}
    \item{\code{SentLength}}{a numeric vector}
    \item{\code{Sessn}}{a numeric vector}
    \item{\code{Set}}{a numeric vector}
    \item{\code{Orgsent}}{a numeric vector}
    \item{\code{TooFast}}{a numeric vector}
    \item{\code{TooFew}}{a numeric vector}
    \item{\code{TooMany}}{a numeric vector}
    \item{\code{TooSlow}}{a numeric vector}
    \item{\code{Wait1.OffsetTime}}{a numeric vector}
    \item{\code{Wait1.OnsetTime}}{a numeric vector}
  }
}
\source{
Johns, C.L., Matsuki, K., & Van Dyke, J. A. (2015). Poor readers' retrieval mechanism: Efficient access is not dependent on reading skill. Manuscript submitted for publication.
}

\examples{
data(Auditory_demo)

}
\keyword{datasets}
