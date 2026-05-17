/*
   @brief Constructs a 
*/
void EventGenerator::SetEventID(int aEventIter)  
{
  ostringstream oss;
  Crc32 tEventID;

  oss << sTimeStamp.Data() << "Event: " << aEventIter;
  tEventID.Update(oss.str().data(), oss.str().length());
  tEventID.Finish(); 

  mEventID = tEventID.GetValue();
}

/*
   @brief Returns the EventID
*/
unsigned int EventGenerator::GetEventID() const
{
  return mEventID;
}
