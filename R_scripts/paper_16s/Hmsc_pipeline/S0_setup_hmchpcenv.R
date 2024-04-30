system_python  <-  "python"
system2(system_python, "--version")
system2(system_python, "-m venv ./R_scripts/paper_16s/Hmsc_pipeline/hmsc-venv")
python <- file.path("./R_scripts/paper_16s/Hmsc_pipeline/hmsc-venv", "Scripts", "python") 
system2(python, "--version")        

package_path <- file.path("./R_scripts/paper_16s/Hmsc_pipeline/require_setup/")
system2(python, "-m pip install --upgrade pip")
system2(python, paste("-m pip install", shQuote(package_path)))


Sys.setenv(TF_CPP_MIN_LOG_LEVEL=3)  # reduce debug output from tensorflow
system2(python, "-c \"import tensorflow as tf; print(tf.constant(1))\"")
system2(python, "-c \"import hmsc\"")
