NETPHOREST_INPUT := /mnt/data/NetPhorest/v2.1/NetPhorest_INPUT


netphorest: insr_code.h nn_code.h nn_data.h pssm_code.h pssm_data.h netphorest.c
	cc -O3 -o netphorest netphorest.c -lm

release: netphorest .public
	ln -f motifs.c netphorest.c
	zip binary_human_release.zip netphorest.c *.h README.txt

clean:
	rm -f binary_human_release.zip
	rm -f netphorest.c

# Deep clean that regenerates all the *.h file
clean-all:
	rm -f insr_code.h nn_code.h nn_data.h pssm_code.h pssm_data.h netphorest.c binary_human_release.zip .public .private

netphorest.c: motifs.c
	ln -f motifs.c netphorest.c

nn:
	./nn2h.pl

insr:
	./insr2h.pl
      
pssm:
	./pssm2h.pl

insr_code.h: insr2h.pl
	./insr2h.pl

nn_code.h nn_data.h: nn2h.pl
	./nn2h.pl

pssm_code.h pssm_data.h: pssm2h.pl
	./pssm2h.pl public $(NETPHOREST_INPUT)/PSSMs

public: pssm_code.h pssm_data.h .public pssm2h.pl

private: pssm_code.h pssm_data.h .private pssm2h.pl

.public: pssm2h.pl
	rm -f .private
	./pssm2h.pl public $(NETPHOREST_INPUT)/PSSMs
	touch .public

.private: pssm2h.pl
	rm -f .public
	./pssm2h.pl private $(NETPHOREST_INPUT)/PSSMs
	touch .private
