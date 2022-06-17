# SOSP_VASO
Stack of Spirals SS-SI-VASO

To perform the SOSP reconstruction:

  1. Clone this repository
  2. Pull docker Julia docker with MRIReco.jl installed using the following command: docker pull monreal93/julia_mri_recon:latest
  3. Start the cointainer with the following command:
     sudo docker run -t -d -P -v /mnt/5T3/Alejandro:/usr/share/Alejandro --name julia_mri_recon monreal93/julia_mri_recon 
     replace "/mnt/5T3/Alejandro" with local folder where data downloaded from Dropbox is located
     replace "/usr/share/Alejandro" with "/usr/share/your_name"
  4. Download the example data from: https://drive.google.com/drive/folders/1kmRqoxQpbzXKYCHgz3m1JStCdjYlmXfE?usp=sharing
  
  
