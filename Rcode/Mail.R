

### R MAILING ###

mail <- function(address, subject, message) {
system(paste("echo '", message,
"' | mail -s '", subject,
"' ", address, sep=""))
}

#We can now send mail from R via, for example,
mail("andrewr", "Test", "Hello world!")

alert <- function() {
mail("andrewr", "Stopped", "Problem")
browser()
}

options(error = alert)

# Alternatively, you can use mutt with exactly the
# same command line structure as mail, if mutt is installed
# on your system. An advantage of mutt is
# that it uses the MIME protocol for binary attachments.
# That would enable you to, say, attach to your email
# the PDF that your script has just created with Sweave
# and pdflatex, or the cvs file that your script creates
# or even the relevant objects saved from the session,
# neatly packaged and compressed in a *.RData object.
# Different flavours of Unix mail are available. You
# can find out more about yours by man mail.

