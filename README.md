# SOSP_VASO
Stack of Spirals SS-SI-VASO

Prerequisits:

Install Docker Engine (https://docs.docker.com/engine/install/ubuntu/)

To perform the SOSP reconstruction:

  1. Clone this repository:
    ```bash
      git clone https://github.com/monreal93/SOSP_VASO.git your_path/sosp_vaso/
    ```
  3. Download the example data from: https://drive.google.com/drive/folders/1L4Dyr4P-co6K44CBIfmMfJGe2vqBTgOR?usp=sharing and save data folder into your_path/sosp_vaso/
  4. Pull docker Julia docker with MRIReco.jl installed using the following command: docker pull monreal93/julia_mri_recon:latest
  5. Start the cointainer with the following command:
     ```
     sudo docker run -t -d -P -v your_path/sosp_vaso:/usr/share/sosp_vaso --name julia_mri_recon monreal93/julia_mri_recon 
     # replace "/mnt/5T3/Alejandro/sosp_vaso" with local folder where data downloaded from Dropbox is located
     ```
  5. Create julia startup file:
     ```bash
     sudo docker exec -it julia_mri_recon bash
     apt-get update
     apt-get install nano
     cd ~
     mkdir .julia/config/
     nano .julia/config/startup.jl
     # In nano editor
     cd("/usr/share/sosp_vaso/recon")  
     using Pkg
     if isfile("Project.toml") && isfile("Manifest.toml") 
      Pkg.activate(".") 
     end
     #
     exit
     ```
  6. Enter julia and instal required packages:
  ```bash
  sudo docker exec -it julia_mri_recon julia -t auto
  Pkg.instantiate()
  ```
  
 #  Perform the reconstruction
 1. you can open and modify the file /sosp_vaso/recon/reconstructions.jl and modify the reconstruction parameters (lines 11-23)
 2. Inside the container run:
 ```julia
 include('reconstructions.jl')
 ```
    
  
  
