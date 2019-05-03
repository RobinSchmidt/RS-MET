
void ImageUpdater::addImageUpdateListener(ImageUpdateListener* listenerToAdd)
{
  imageUpdateListeners.addIfNotAlreadyThere(listenerToAdd);
}

void ImageUpdater::removeImageUpdateListener(ImageUpdateListener* listenerToRemove)
{
  imageUpdateListeners.removeFirstMatchingValue(listenerToRemove);
}

void ImageUpdater::sendImageUpdateNotification(Image* image)
{
  for(int i = 0; i < imageUpdateListeners.size(); i++)
    imageUpdateListeners[i]->imageWasUpdated(image);
}