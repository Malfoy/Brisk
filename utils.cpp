

//START LOW LEVEL FUNCTIONS



string getLineFasta(zstr::ifstream* in) {
	string line, result;
	getline(*in, line);
	char c = static_cast<char>(in->peek());
	while (c != '>' and c != EOF) {
		getline(*in, line);
		result += line;
		c = static_cast<char>(in->peek());
	}
	return result;
}



uint32_t revhash(uint32_t x) {
	x = ((x >> 16) ^ x) * 0x2c1b3c6d;
	x = ((x >> 16) ^ x) * 0x297a2d39;
	x = ((x >> 16) ^ x);
	return x;
}



uint32_t unrevhash(uint32_t x) {
	x = ((x >> 16) ^ x) * 0x0cf0b109; // PowerMod[0x297a2d39, -1, 2^32]
	x = ((x >> 16) ^ x) * 0x64ea2d65;
	x = ((x >> 16) ^ x);
	return x;
}



uint64_t revhash(uint64_t x) {
	x = ((x >> 32) ^ x) * 0xD6E8FEB86659FD93;
	x = ((x >> 32) ^ x) * 0xD6E8FEB86659FD93;
	x = ((x >> 32) ^ x);
	return x;
}



uint64_t unrevhash(uint64_t x) {
	x = ((x >> 32) ^ x) * 0xCFEE444D8B59A89B;
	x = ((x >> 32) ^ x) * 0xCFEE444D8B59A89B;
	x = ((x >> 32) ^ x);
	return x;
}



uint64_t nuc2int(char c) {
	return (c / 2) % 4;
}



uint64_t nuc2intrc(char c) {
	return ((c / 2) % 4) ^ 2;
}



void updateK(uint64_t& min, char nuc) {
	min <<= 2;
	min += nuc2int(nuc);
	min %= offsetUpdateAnchor;
}



 void add_nuc_superkmer(SKC& min, char nuc) {
	min.sk <<= 2;
	min.sk += nuc2int(nuc);
	min.counts[min.size++] = 0;
}



 void updateM(uint64_t& min, char nuc) {
	min <<= 2;
	min += nuc2int(nuc);
	min %= offsetUpdateAnchorMin;
}



void updateRCK(uint64_t& min, char nuc) {
	min >>= 2;
	min += (nuc2intrc(nuc) << (2 * k - 2));
}



void updateRCM(uint64_t& min, char nuc) {
	min >>= 2;
	min += (nuc2intrc(nuc) << (2 * minimizer_size - 2));
}



char revCompChar(char c) {
	switch (c) {
		case 'A':
			return 'T';
		case 'C':
			return 'G';
		case 'G':
			return 'C';
	}
	return 'A';
}



string revComp(const string& s) {
	string rc(s.size(), 0);
	for (int i((int)s.length() - 1); i >= 0; i--) {
		rc[s.size() - 1 - i] = revCompChar(s[i]);
	}
	return rc;
}



string getCanonical(const string& str) {
	return (min(str, revComp(str)));
}



string intToString(uint64_t n) {
	if (n < 1000) {
		return to_string(n);
	}
	string end(to_string(n % 1000));
	if (end.size() == 3) {
		return intToString(n / 1000) + "," + end;
	}
	if (end.size() == 2) {
		return intToString(n / 1000) + ",0" + end;
	}
	return intToString(n / 1000) + ",00" + end;
}



//START HIGH LEVEL FUNCTIONS
