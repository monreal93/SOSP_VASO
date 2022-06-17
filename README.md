# SOSP_VASO
Stack of Spirals SS-SI-VASO

Prerequisits:

Install Docker Engine (https://docs.docker.com/engine/install/ubuntu/)

To perform the SOSP reconstruction:

  1. Clone this repository
  2. Download the example data from: https://drive.google.com/drive/folders/1kmRqoxQpbzXKYCHgz3m1JStCdjYlmXfE?usp=sharing
  3. Pull docker Julia docker with MRIReco.jl installed using the following command: docker pull monreal93/julia_mri_recon:latest
  4. Start the cointainer with the following command:
     ```
     sudo docker run -t -d -P -v /mnt/5T3/Alejandro/sosp_vaso:/usr/share/sosp_vaso --name julia_mri_recon monreal93/julia_mri_recon 
     # replace "/mnt/5T3/Alejandro/sosp_vaso" with local folder where data downloaded from Dropbox is located
     ```
  5. Create julia startup file:
     ``` bash
     sudo docker exec -it julia_mri_recon bash
     apt-get update
     apt-get install nano
     cd ~
     mkdir .julia/config/
     nano .julia/config/startup.jl
     # In editor
     cd("/usr/share/sosp_vaso/recon")  
     using Pkg
     if isfile("Project.toml") && isfile("Manifest.toml") 
      Pkg.activate(".") 
     end 
     ```
  6. Enter julia:
    ```
    julia -t auto
    ```
  7. Install required packages: 
    ```
    Pkg.instantiate()
    ```
    
  
  
