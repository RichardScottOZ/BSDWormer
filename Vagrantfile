# -*- mode: ruby -*-
# vi: set ft=ruby :

$provision_script = <<'EOF'
echo "shell provisioning"
PG_VERSION=9.3
PROVISIONED_ON=/etc/vm_provision_on_timestamp
sudo apt-get update
sudo apt-get -y dist-upgrade
sudo apt-get install -y virtualbox-guest-dkms python3-scipy ipython3-notebook python3-matplotlib git python3-pip python3-gdal gdal-bin python3-networkx atop python3-sqlalchemy postgresql-9.3-postgis-2.1 python3-psycopg2
sudo -u postgres createuser vagrant
pip3 install future
pip3 install geoalchemy2
pip3 install git+https://fghorow@bitbucket.org/fghorow/pyvtk-py3-port.git

PG_CONF="/etc/postgresql/$PG_VERSION/main/postgresql.conf"
PG_HBA="/etc/postgresql/$PG_VERSION/main/pg_hba.conf"
PG_DIR="/var/lib/postgresql/$PG_VERSION/main"

# Edit postgresql.conf to change listen address to '*':
sed -i "s/#listen_addresses = 'localhost'/listen_addresses = '*'/" "$PG_CONF"

# Append to pg_hba.conf to add password auth:
echo "host    all             all             all                     md5" >> "$PG_HBA"

# Explicitly set default client_encoding
echo "client_encoding = utf8" >> "$PG_CONF"

# Restart so that all new config is loaded:
service postgresql restart

# Edit the following to change the name of the database user that will be created:
APP_DB_USER=wormuser
APP_DB_PASS=wormuserpass

# Edit the following to change the name of the database that is created (defaults to the user name)
APP_DB_NAME=$APP_DB_USER

cat << EOF1 | su - postgres -c psql
-- Create the database user:
CREATE USER $APP_DB_USER WITH PASSWORD '$APP_DB_PASS';

-- Create the database:
CREATE DATABASE $APP_DB_NAME WITH OWNER=$APP_DB_USER
                                  LC_COLLATE='en_US.utf8'
                                  LC_CTYPE='en_US.utf8'
                                  ENCODING='UTF8'
                                  TEMPLATE=template0;
EOF1

# This enables the PostGIS stuff inside the right database
sudo -u postgres psql -c "CREATE EXTENSION postgis; CREATE EXTENSION postgis_topology;" $APP_DB_NAME

# Tag the provision time:
date > "$PROVISIONED_ON"



echo "Successfully created PostgreSQL dev virtual machine."
echo ""

EOF

# All Vagrant configuration is done below. The "2" in Vagrant.configure
# configures the configuration version (we support older styles for
# backwards compatibility). Please don't change it unless you know what
# you're doing.
Vagrant.configure(2) do |config|
  # The most common configuration options are documented and commented below.
  # For a complete reference, please see the online documentation at
  # https://docs.vagrantup.com.

  # Every Vagrant development environment requires a box. You can search for
  # boxes at https://atlas.hashicorp.com/search.
  config.vm.box = "ubuntu/trusty64"

  # Disable automatic box update checking. If you disable this, then
  # boxes will only be checked for updates when the user runs
  # `vagrant box outdated`. This is not recommended.
  # config.vm.box_check_update = false

  # Create a forwarded port mapping which allows access to a specific port
  # within the machine from a port on the host machine. In the example below,
  # accessing "localhost:8080" will access port 80 on the guest machine.
  config.vm.network "forwarded_port", guest: 8888, host: 8888, auto_correct: true
 
  config.vm.network "forwarded_port", guest: 5432, host: 15432, auto_correct: true
 
  # Create a private network, which allows host-only access to the machine
  # using a specific IP.
  config.vm.network "private_network", ip: "192.168.33.10"

  # Create a public network, which generally matched to bridged network.
  # Bridged networks make the machine appear as another physical device on
  # your network.
  # config.vm.network "public_network"

  # Share an additional folder to the guest VM. The first argument is
  # the path on the host to the actual folder. The second argument is
  # the path on the guest to mount the folder. And the optional third
  # argument is a set of non-required options.
  config.vm.synced_folder "./", "/vagrant"

  # Provider-specific configuration so you can fine-tune various
  # backing providers for Vagrant. These expose provider-specific options.
  # Example for VirtualBox:
  #
  config.vm.provider "virtualbox" do |vb|
    # Display the VirtualBox GUI when booting the machine
    # vb.gui = true
    
    # Customize the amount of memory on the VM:
    vb.memory = "4096"
  end
  #
  # View the documentation for the provider you are using for more
  # information on available options.

  # Define a Vagrant Push strategy for pushing to Atlas. Other push strategies
  # such as FTP and Heroku are also available. See the documentation at
  # https://docs.vagrantup.com/v2/push/atlas.html for more information.
  # config.push.define "atlas" do |push|
  #   push.app = "YOUR_ATLAS_USERNAME/YOUR_APPLICATION_NAME"
  # end
  
  
  # we will try to autodetect this path. 
  # However, if we cannot or you have a special one you may pass it like:
  # config.vbguest.iso_path = "#{ENV['HOME']}/Downloads/VBoxGuestAdditions.iso"
  # or
  # config.vbguest.iso_path = "http://company.server/VirtualBox/%{version}/VBoxGuestAdditions.iso"

  # set auto_update to false, if you do NOT want to check the correct 
  # additions version when booting this machine
  config.vbguest.auto_update = false

  # do NOT download the iso file from a webserver
  config.vbguest.no_remote = true

  # Enable provisioning with a shell script. Additional provisioners such as
  # Puppet, Chef, Ansible, Salt, and Docker are also available. Please see the
  # documentation for more information about their specific syntax and use.
  config.vm.provision :shell, :inline => $provision_script
  
  # Fire off the Ipython notebook server
  config.vm.provision "shell", run: "always", inline: <<-SHELL
    ipython3 notebook --notebook-dir=/vagrant/src --no-browser --ip=0.0.0.0 &
  SHELL

end
