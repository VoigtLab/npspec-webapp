# metabs-webapp
Webapp to display predicted spectral data for metabolites

For simple deployment: `streamlit run Home.py`


To deploy a website that can be accessed with https :
- Obtain SSL certificates for your server using [certbot](https://certbot.eff.org/instructions?ws=nginx&os=pip): `certbot --nginx -d <SERVER NAME>`
- Copy certificates to this project directory, e.g. (need both):
	- cp /etc/letsencrypt/live/<SERVER NAME>/fullchain.pem ./nginx/certs/<SERVER NAME>.crt`
	- cp /etc/letsencrypt/live/<SERVER NAME>/privkey.pem nginx/certs/<SERVER NAME>.key
- Modify nginx config file:
	- Copy template file `cp default_nginx_template.conf nginx/conf.d/default.conf`
	- Edit the template file to replace <SERVER-NAME> with your server name
- Deploy with Docker:
	- `docker compose up -d`
