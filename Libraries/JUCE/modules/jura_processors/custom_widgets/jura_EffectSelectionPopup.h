#ifndef jura_EffectSelectionPopup_h
#define jura_EffectSelectionPopup_h

/**  This is Popup menu to select an effect-algorithm for Quadrifex. 

\todo: rename to QuadrifexEffectMenu and maybe move into Quadrifex.h/cpp - it's not used 
elsewhere. */

class EffectSelectionPopup : public RPopUpMenu
{

public:

  //---------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  EffectSelectionPopup(Component *componentToAttachTo);

protected:

  /*
  // the sub-menus (obsolete - use a single menu with treeview now):
  RPopUpMenuOld stereoToolsMenu, filterMenu, delayMenu, distortionMenu, modulationMenu,
    dynamicsMenu, fatteningMenu, spectralMenu, miscMenu, generatorsMenu;
  // Delays, Dynamics, Spectral, Modulation,....
  */

  /** Opens the PopupMenu that appears on right clicks. */
  //void openRightClickPopupMenu();

  ///void createRightClickPopupMenu(PopupMenu*& mainMenu, PopupMenu*& defaultValueSubMenu,
  //  int& defaultValueIndicesMin, int& defaultValueIndicesMax);

  /** Handles the result of opening the right-click popoup menu. */
  //void handleRightClickPopupMenuResult(int result, int defaultValueIndicesMin, int defaultValueIndicesMax);

  juce_UseDebuggingNewOperator;
};

#endif  
