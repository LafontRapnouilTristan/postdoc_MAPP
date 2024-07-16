system_python <- "python"
system2(system_python, "--version")

system2(system_python, "-m venv hmsc-venv")

python <- file.path(getwd(), "hmsc-venv", "Scripts", "python")  # for Windows
system2(python, "--version")


python <- file.path(getwd(), "hmsc-venv", "Scripts", "python")  # for Windows
system2(python, "--version")

package_path <- file.path(getwd(), "..", "..")
system2(python, "-m pip install --upgrade pip")
system2(python, paste("-m pip install", shQuote(package_path)))

python <- file.path(getwd(), "hmsc-venv", "Scripts", "python")  # hmsc-venv for Windows

# reticulate::py_install("tensorflow","./hmsc-venv")

Sys.setenv(TF_CPP_MIN_LOG_LEVEL=3)  # reduce debug output from tensorflow
system2(python, "-c \"import tensorflow as tf; print(tf.constant(1))\"")
system2(python, "-c \"import hmsc\"")


library(devtools)
install_github("hmsc-r/HMSC")

