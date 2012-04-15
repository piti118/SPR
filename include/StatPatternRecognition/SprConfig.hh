//$Id: SprConfig.hh,v 1.1 2008-05-08 19:57:43 narsky Exp $
/* 
* Author: Piti Ongmongkolkul
* this class is intended for reading config file of the form
* xxx=yyy
* zzz=aaa
*/
#ifndef SprConfig_HH
#define SprConfig_HH

#include <map>
#include <string>
#include <iostream>

class SprConfig
{
private:
  std::map<std::string,std::string> config_;

  static const char commentSignal_;
  static const char keyValueSeparator_;

  static void trim(std::string& s, const std::string drop=" ");
  //  static std::string ltrim(std::string& s,const std::string drop=" ");

public:
  ~SprConfig() {}

  SprConfig(const std::string& filename);
  SprConfig(const SprConfig& configObject);  
  
  // read config from file
  bool read(const std::string& filename);

  //accessor with conversion
  std::string getStringValue(const std::string& key, 
			     const std::string defaultValue=" ") const;
  int getIntValue(const std::string& key, const int defaultValue=0) const;
  double getDoubleValue(const std::string& key, 
			const double defaultValue=0.0) const;
  bool getBoolValue(const std::string& key, 
		    const bool defaultValue=false) const;  
  bool keyExists(const std::string& key) const;
  void print(std::ostream& os) const;
};

#endif
