/*
  ==============================================================================

   This file is part of the JUCE library.
   Copyright (c) 2017 - ROLI Ltd.

   JUCE is an open source library subject to commercial or open-source
   licensing.

   By using JUCE, you agree to the terms of both the JUCE 5 End-User License
   Agreement and JUCE 5 Privacy Policy (both updated and effective as of the
   27th April 2017).

   End User License Agreement: www.juce.com/juce-5-licence
   Privacy Policy: www.juce.com/juce-5-privacy-policy

   Or: You may also use this code under the terms of the GPL v3 (see
   www.gnu.org/licenses).

   JUCE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY, AND ALL WARRANTIES, WHETHER
   EXPRESSED OR IMPLIED, INCLUDING MERCHANTABILITY AND FITNESS FOR PURPOSE, ARE
   DISCLAIMED.

  ==============================================================================
*/

namespace juce
{

#define JNI_CLASS_MEMBERS(METHOD, STATICMETHOD, FIELD, STATICFIELD, CALLBACK) \
    METHOD (isBillingSupported,      "isBillingSupported",      "(ILjava/lang/String;Ljava/lang/String;)I") \
    METHOD (getSkuDetails,           "getSkuDetails",           "(ILjava/lang/String;Ljava/lang/String;Landroid/os/Bundle;)Landroid/os/Bundle;") \
    METHOD (getBuyIntent,            "getBuyIntent",            "(ILjava/lang/String;Ljava/lang/String;Ljava/lang/String;Ljava/lang/String;)Landroid/os/Bundle;") \
    METHOD (getBuyIntentExtraParams, "getBuyIntentExtraParams", "(ILjava/lang/String;Ljava/lang/String;Ljava/lang/String;Ljava/lang/String;Landroid/os/Bundle;)Landroid/os/Bundle;") \
    METHOD (getPurchases,            "getPurchases",            "(ILjava/lang/String;Ljava/lang/String;Ljava/lang/String;)Landroid/os/Bundle;") \
    METHOD (consumePurchase,         "consumePurchase",         "(ILjava/lang/String;Ljava/lang/String;)I") \
    METHOD (getPurchaseHistory,      "getPurchaseHistory",      "(ILjava/lang/String;Ljava/lang/String;Ljava/lang/String;Landroid/os/Bundle;)Landroid/os/Bundle;")

DECLARE_JNI_CLASS (IInAppBillingService, "com/android/vending/billing/IInAppBillingService")
#undef JNI_CLASS_MEMBERS

#define JNI_CLASS_MEMBERS(METHOD, STATICMETHOD, FIELD, STATICFIELD, CALLBACK) \
    STATICMETHOD (asInterface,      "asInterface",      "(Landroid/os/IBinder;)Lcom/android/vending/billing/IInAppBillingService;") \

DECLARE_JNI_CLASS (IInAppBillingServiceStub, "com/android/vending/billing/IInAppBillingService$Stub")
#undef JNI_CLASS_MEMBERS

//==============================================================================
struct ServiceConnection  : public AndroidInterfaceImplementer
{
    virtual void onServiceConnected    (jobject component, jobject iBinder) = 0;
    virtual void onServiceDisconnected (jobject component) = 0;

    jobject invoke (jobject proxy, jobject method, jobjectArray args) override
    {
        auto* env = getEnv();
        auto methodName = juceString ((jstring) env->CallObjectMethod (method, JavaMethod.getName));

        if (methodName == "onServiceConnected")
        {
            onServiceConnected (env->GetObjectArrayElement (args, 0),
                                env->GetObjectArrayElement (args, 1));
            return nullptr;
        }

        if (methodName == "onServiceDisconnected")
        {
            onServiceDisconnected (env->GetObjectArrayElement (args, 0));
            return nullptr;
        }

        return AndroidInterfaceImplementer::invoke (proxy, method, args);
    }
};

//==============================================================================
struct InAppPurchases::Pimpl    : private AsyncUpdater,
                                  private ServiceConnection
{
    Pimpl (InAppPurchases& parent)  : owner (parent)
    {
        auto* env = getEnv();
        auto intent = env->NewObject (AndroidIntent, AndroidIntent.constructWithString,
                                      javaString ("com.android.vending.billing.InAppBillingService.BIND").get());
        env->CallObjectMethod (intent, AndroidIntent.setPackage, javaString ("com.android.vending").get());

        serviceConnection = GlobalRef (CreateJavaInterface (this, "android/content/ServiceConnection"));

        env->CallBooleanMethod (getCurrentActivity().get(), AndroidContext.bindService, intent,
                                serviceConnection.get(), 1 /*BIND_AUTO_CREATE*/);

        if (threadPool == nullptr)
            threadPool.reset (new ThreadPool (1));
    }

    ~Pimpl()
    {
        threadPool = nullptr;

        if (serviceConnection != nullptr)
        {
            getEnv()->CallVoidMethod (getCurrentActivity().get(), AndroidContext.unbindService, serviceConnection.get());
            serviceConnection.clear();
        }
    }

    //==============================================================================
    bool isInAppPurchasesSupported()       { return isInAppPurchasesSupported (inAppBillingService); }

    void getProductsInformation (const StringArray& productIdentifiers)
    {
        auto callback = [this](const Array<InAppPurchases::Product>& products)
        {
            const ScopedLock lock (getProductsInformationJobResultsLock);
            getProductsInformationJobResults.insert (0, products);
            triggerAsyncUpdate();
        };

        threadPool->addJob (new GetProductsInformationJob (*this, getPackageName(),
                                                           productIdentifiers, callback), true);
    }

    void purchaseProduct (const String& productIdentifier, bool isSubscription,
                          const StringArray& subscriptionIdentifiers, bool creditForUnusedSubscription)
    {
        // Upgrading/downgrading only makes sense for subscriptions!
        jassert (subscriptionIdentifiers.isEmpty() || isSubscription);

        auto buyIntentBundle = getBuyIntentBundle (productIdentifier, isSubscription,
                                                   subscriptionIdentifiers, creditForUnusedSubscription);
        auto* env = getEnv();

        auto responseCodeString = javaString ("RESPONSE_CODE");
        auto responseCode = env->CallIntMethod (buyIntentBundle.get(), AndroidBundle.getInt, responseCodeString.get());

        if (responseCode == 0)
        {
            auto buyIntentString = javaString ("BUY_INTENT");
            auto pendingIntent   = LocalRef<jobject> (env->CallObjectMethod (buyIntentBundle.get(), AndroidBundle.getParcelable, buyIntentString.get()));

            auto  requestCode = 1001;
            auto intentSender    = LocalRef<jobject> (env->CallObjectMethod (pendingIntent.get(), AndroidPendingIntent.getIntentSender));
            auto fillInIntent    = LocalRef<jobject> (env->NewObject (AndroidIntent, AndroidIntent.constructor));
            auto flagsMask       = LocalRef<jobject> (env->CallStaticObjectMethod (JavaInteger, JavaInteger.valueOf, 0));
            auto flagsValues     = LocalRef<jobject> (env->CallStaticObjectMethod (JavaInteger, JavaInteger.valueOf, 0));
            auto extraFlags      = LocalRef<jobject> (env->CallStaticObjectMethod (JavaInteger, JavaInteger.valueOf, 0));

            env->CallVoidMethod (getCurrentActivity().get(), AndroidActivity.startIntentSenderForResult, intentSender.get(), requestCode,
                                 fillInIntent.get(), flagsMask.get(), flagsValues.get(), extraFlags.get());
        }
        else if (responseCode == 7)
        {
            // Item already bought.
            notifyAboutPurchaseResult ({ {}, productIdentifier, juceString (getPackageName()), {}, {} }, true, statusCodeToUserString (responseCode));
        }
    }

    void restoreProductsBoughtList (bool, const juce::String&)
    {
        auto callback = [this](const GetProductsBoughtJob::Result& r)
        {
            const ScopedLock lock (getProductsBoughtJobResultsLock);
            getProductsBoughtJobResults.insert (0, r);
            triggerAsyncUpdate();
        };

        threadPool->addJob (new GetProductsBoughtJob (*this,
                                                      getPackageName(), callback), true);
    }

    void consumePurchase (const String& productIdentifier, const String& purchaseToken)
    {
        auto callback = [this](const ConsumePurchaseJob::Result& r)
        {
            const ScopedLock lock (consumePurchaseJobResultsLock);
            consumePurchaseJobResults.insert (0, r);
            triggerAsyncUpdate();
        };

        threadPool->addJob (new ConsumePurchaseJob (*this, getPackageName(), productIdentifier,
                                                    purchaseToken, callback), true);
    }

    //==============================================================================
    void startDownloads  (const Array<Download*>& downloads)
    {
        // Not available on this platform.
        ignoreUnused (downloads);
        jassertfalse;
    }

    void pauseDownloads  (const Array<Download*>& downloads)
    {
        // Not available on this platform.
        ignoreUnused (downloads);
        jassertfalse;
    }

    void resumeDownloads  (const Array<Download*>& downloads)
    {
        // Not available on this platform.
        ignoreUnused (downloads);
        jassertfalse;
    }

    void cancelDownloads  (const Array<Download*>& downloads)
    {
        // Not available on this platform.
        ignoreUnused (downloads);
        jassertfalse;
    }

    //==============================================================================
    LocalRef<jobject> getBuyIntentBundle (const String& productIdentifier, bool isSubscription,
                                          const StringArray& subscriptionIdentifiers, bool creditForUnusedSubscription)
    {
        auto* env = getEnv();

        auto skuString         = javaString (productIdentifier);
        auto productTypeString = javaString (isSubscription ? "subs" : "inapp");
        auto devString         = javaString ("");

        if (subscriptionIdentifiers.isEmpty())
            return LocalRef<jobject> (inAppBillingService.callObjectMethod (IInAppBillingService.getBuyIntent, 3,
                                                                            getPackageName().get(), skuString.get(),
                                                                            productTypeString.get(), devString.get()));

        auto skuList = LocalRef<jobject> (env->NewObject (JavaArrayList, JavaArrayList.constructor,
                                                          (int) subscriptionIdentifiers.size()));

        if (skuList.get() == 0)
        {
            jassertfalse;
            return LocalRef<jobject> (0);
        }

        for (const auto& identifier : subscriptionIdentifiers)
            env->CallBooleanMethod (skuList.get(), JavaArrayList.add, javaString (identifier).get());

        auto extraParams = LocalRef<jobject> (env->NewObject (AndroidBundle, AndroidBundle.constructor));

        if (extraParams.get() == 0)
        {
            jassertfalse;
            return LocalRef<jobject> (0);
        }

        auto skusToReplaceString        = javaString ("skusToReplace");
        auto replaceSkusProrationString = javaString ("replaceSkusProration");

        env->CallVoidMethod (extraParams.get(), AndroidBundle.putStringArrayList, skusToReplaceString.get(), skuList.get());
        env->CallVoidMethod (extraParams.get(), AndroidBundle.putBoolean, replaceSkusProrationString.get(), creditForUnusedSubscription);

        return LocalRef<jobject> (inAppBillingService.callObjectMethod (IInAppBillingService.getBuyIntentExtraParams, 6,
                                                                        getPackageName().get(), skuString.get(),
                                                                        productTypeString.get(), devString.get(),
                                                                        extraParams.get()));
    }

    //==============================================================================
    void notifyAboutPurchaseResult (const InAppPurchases::Purchase& purchase, bool success, const String& statusDescription)
    {
        owner.listeners.call ([&] (Listener& l) { l.productPurchaseFinished ({ purchase, {} }, success, statusDescription); });
    }

    //==============================================================================
    bool checkIsReady()
    {
        // It may take a few seconds for the in-app purchase service to connect
        for (auto retries = 0; retries < 10 && inAppBillingService.get() == 0; ++retries)
            Thread::sleep (500);

        return (inAppBillingService.get() != 0);
    }

    static bool isInAppPurchasesSupported (jobject iapService)
    {
        if (iapService != nullptr)
        {
            auto* env = getEnv();

            auto inAppString = javaString ("inapp");
            auto subsString  = javaString ("subs");

            if (env->CallIntMethod (iapService, IInAppBillingService.isBillingSupported, 3,
                                    getPackageName().get(), inAppString.get()) != 0)
                return false;

            if (env->CallIntMethod (iapService, IInAppBillingService.isBillingSupported, 3,
                                    getPackageName().get(), subsString.get()) != 0)
                return false;

            return true;
        }

        // Connecting to the in-app purchase server failed! This could have multiple reasons:
        // 1) Your phone/emulator must support the google play store
        // 2) Your phone must be logged into the google play store and be able to receive updates
        // 3) It can take a few seconds after instantiation of the InAppPurchase class for
        //    in-app purchases to be avaialable on Android.
        return false;
    }

    //==============================================================================
    void onServiceConnected (jobject, jobject iBinder) override
    {
        auto* env = getEnv();

        LocalRef<jobject> iapService (env->CallStaticObjectMethod (IInAppBillingServiceStub,
                                                                   IInAppBillingServiceStub.asInterface,
                                                                   iBinder));

        if (isInAppPurchasesSupported (iapService))
            inAppBillingService = GlobalRef (iapService);

        // If you hit this assert, then in-app purchases is not available on your device,
        // most likely due to too old version of Google Play API (hint: update Google Play on the device).
        jassert (isInAppPurchasesSupported());
    }

    void onServiceDisconnected (jobject) override
    {
        inAppBillingService.clear();
    }

    //==============================================================================
    static LocalRef<jstring> getPackageName()
    {
        return LocalRef<jstring> ((jstring) (getEnv()->CallObjectMethod (getAppContext().get(), AndroidContext.getPackageName)));
    }

    //==============================================================================
    struct GetProductsInformationJob  : public ThreadPoolJob
    {
        using Callback = std::function<void(const Array<InAppPurchases::Product>&)>;

        GetProductsInformationJob (Pimpl& parent,
                                   const LocalRef<jstring>& packageNameToUse,
                                   const StringArray& productIdentifiersToUse,
                                   const Callback& callbackToUse)
            : ThreadPoolJob ("GetProductsInformationJob"),
              owner (parent),
              packageName (LocalRef<jobject> (packageNameToUse.get())),
              productIdentifiers (productIdentifiersToUse),
              callback (callbackToUse)
        {}

        ThreadPoolJob::JobStatus runJob() override
        {
            jassert (callback);

            if (owner.checkIsReady())
            {
                // Google's Billing API limitation
                auto maxQuerySize = 20;
                auto pi = 0;

                Array<InAppPurchases::Product> results;
                StringArray identifiersToUse;

                for (auto i = 0; i < productIdentifiers.size(); ++i)
                {
                    identifiersToUse.add (productIdentifiers[i].toLowerCase());
                    ++pi;

                    if (pi == maxQuerySize || i == productIdentifiers.size() - 1)
                    {
                        auto inAppProducts = processRetrievedProducts (queryProductsInformationFromService (identifiersToUse, "inapp"));
                        auto subsProducts  = processRetrievedProducts (queryProductsInformationFromService (identifiersToUse, "subs"));

                        results.addArray (inAppProducts);
                        results.addArray (subsProducts);
                        identifiersToUse.clear();
                        pi = 0;
                    }
                }

                if (callback)
                    callback (results);
            }
            else
            {
                if (callback)
                    callback ({});
            }

            return jobHasFinished;
        }

    private:
        LocalRef<jobject> queryProductsInformationFromService (const StringArray& productIdentifiersToQuery, const String& productType)
        {
            auto* env = getEnv();

            auto skuList = LocalRef<jobject> (env->NewObject (JavaArrayList, JavaArrayList.constructor, productIdentifiersToQuery.size()));

            if (skuList.get() == 0)
                return LocalRef<jobject> (0);

            for (const auto& pi : productIdentifiersToQuery)
                env->CallBooleanMethod (skuList.get(), JavaArrayList.add, javaString (pi).get());

            auto querySkus = LocalRef<jobject> (env->NewObject (AndroidBundle, AndroidBundle.constructor));

            if (querySkus.get() == 0)
                return LocalRef<jobject> (0);

            auto itemIdListString = javaString ("ITEM_ID_LIST");

            env->CallVoidMethod (querySkus.get(), AndroidBundle.putStringArrayList, itemIdListString.get(), skuList.get());

            auto productTypeString = javaString (productType);

            auto productDetails = LocalRef<jobject> (owner.inAppBillingService.callObjectMethod (IInAppBillingService.getSkuDetails,
                                                                                           3, (jstring) packageName.get(),
                                                                                           productTypeString.get(), querySkus.get()));

            return productDetails;
        }

        Array<InAppPurchases::Product> processRetrievedProducts (LocalRef<jobject> retrievedProducts)
        {
            Array<InAppPurchases::Product> products;

            if (owner.checkIsReady())
            {
                auto* env = getEnv();

                auto responseCodeString = javaString ("RESPONSE_CODE");

                auto responseCode = env->CallIntMethod (retrievedProducts.get(), AndroidBundle.getInt, responseCodeString.get());

                if (responseCode == 0)
                {
                    auto detailsListString = javaString ("DETAILS_LIST");

                    auto responseList = LocalRef<jobject> (env->CallObjectMethod (retrievedProducts.get(), AndroidBundle.getStringArrayList,
                                                                                  detailsListString.get()));

                    if (responseList != 0)
                    {
                        auto iterator = LocalRef<jobject> (env->CallObjectMethod (responseList.get(), JavaArrayList.iterator));

                        if (iterator.get() != 0)
                        {
                            for (;;)
                            {
                                if (! env->CallBooleanMethod (iterator, JavaIterator.hasNext))
                                    break;

                                auto response = juce::LocalRef<jstring> ((jstring)env->CallObjectMethod (iterator, JavaIterator.next));

                                if (response.get() != 0)
                                {
                                    var responseData = JSON::parse (juceString (response.get()));

                                    if (DynamicObject* object = responseData.getDynamicObject())
                                    {
                                        NamedValueSet& props = object->getProperties();

                                        static Identifier productIdIdentifier         ("productId");
                                        static Identifier titleIdentifier             ("title");
                                        static Identifier descriptionIdentifier       ("description");
                                        static Identifier priceIdentifier             ("price");
                                        static Identifier priceCurrencyCodeIdentifier ("price_currency_code");

                                        var productId         = props[productIdIdentifier];
                                        var title             = props[titleIdentifier];
                                        var description       = props[descriptionIdentifier];
                                        var price             = props[priceIdentifier];
                                        var priceCurrencyCode = props[priceCurrencyCodeIdentifier];

                                        products.add ( { productId.toString(),
                                                         title.toString(),
                                                         description.toString(),
                                                         price.toString(),
                                                         priceCurrencyCode.toString() } );
                                    }

                                }
                            }
                        }
                    }
                }
            }

            return products;
        }

        Pimpl& owner;
        GlobalRef packageName;
        const StringArray productIdentifiers;
        Callback callback;
    };

    //==============================================================================
    struct GetProductsBoughtJob  : public ThreadPoolJob
    {
        struct Result
        {
            bool success = false;
            Array<InAppPurchases::Listener::PurchaseInfo> purchases;
            String statusDescription;
        };

        using Callback = std::function<void(const Result&)>;

        GetProductsBoughtJob (Pimpl& parent,
                              const LocalRef<jstring>& packageNameToUse,
                              const Callback& callbackToUse)
            : ThreadPoolJob ("GetProductsBoughtJob"),
              owner (parent),
              packageName (LocalRef<jobject> (packageNameToUse.get())),
              callback (callbackToUse)
        {}

        ThreadPoolJob::JobStatus runJob() override
        {
            jassert (callback);

            if (owner.checkIsReady())
            {
                auto inAppPurchases = getProductsBought ("inapp", 0);
                auto subsPurchases  = getProductsBought ("subs", 0);

                inAppPurchases.addArray (subsPurchases);

                Array<InAppPurchases::Listener::PurchaseInfo> purchases;

                for (const auto& purchase : inAppPurchases)
                    purchases.add ({ purchase, {} });

                if (callback)
                    callback ({true, purchases, "Success"});
            }
            else
            {
                if (callback)
                    callback ({false, {}, "In-App purchases unavailable"});
            }

            return jobHasFinished;
        }

    private:
        Array<InAppPurchases::Purchase> getProductsBought (const String& productType, jstring continuationToken)
        {
            Array<InAppPurchases::Purchase> purchases;
            auto* env = getEnv();

            auto productTypeString = javaString (productType);
            auto ownedItems = LocalRef<jobject> (owner.inAppBillingService.callObjectMethod (IInAppBillingService.getPurchases, 3,
                                                                                       (jstring) packageName.get(), productTypeString.get(),
                                                                                       continuationToken));

            if (ownedItems.get() != 0)
            {
                auto responseCodeString = javaString ("RESPONSE_CODE");
                auto responseCode = env->CallIntMethod (ownedItems.get(), AndroidBundle.getInt, responseCodeString.get());

                if (responseCode == 0)
                {
                    auto itemListString          = javaString ("INAPP_PURCHASE_ITEM_LIST");
                    auto dataListString          = javaString ("INAPP_PURCHASE_DATA_LIST");
                    auto signatureListString     = javaString ("INAPP_DATA_SIGNATURE_LIST");
                    auto continuationTokenString = javaString ("INAPP_CONTINUATION_TOKEN");

                    auto ownedSkus            = LocalRef<jobject> (env->CallObjectMethod (ownedItems.get(), AndroidBundle.getStringArrayList, itemListString.get()));
                    auto purchaseDataList     = LocalRef<jobject> (env->CallObjectMethod (ownedItems.get(), AndroidBundle.getStringArrayList, dataListString.get()));
                    auto signatureList        = LocalRef<jobject> (env->CallObjectMethod (ownedItems.get(), AndroidBundle.getStringArrayList, signatureListString.get()));
                    auto newContinuationToken = LocalRef<jstring> ((jstring) env->CallObjectMethod (ownedItems.get(), AndroidBundle.getString, continuationTokenString.get()));

                    for (auto i = 0; i < env->CallIntMethod (purchaseDataList.get(), JavaArrayList.size); ++i)
                    {
                        auto sku          = juceString ((jstring) (env->CallObjectMethod (ownedSkus.get(),        JavaArrayList.get, i)));
                        auto purchaseData = juceString ((jstring) (env->CallObjectMethod (purchaseDataList.get(), JavaArrayList.get, i)));
                        auto signature    = juceString ((jstring) (env->CallObjectMethod (signatureList.get(),    JavaArrayList.get, i)));

                        var responseData = JSON::parse (purchaseData);

                        if (auto* object = responseData.getDynamicObject())
                        {
                            auto& props = object->getProperties();

                            static const Identifier orderIdIdentifier       ("orderId"),
                                                    packageNameIdentifier   ("packageName"),
                                                    productIdIdentifier     ("productId"),
                                                    purchaseTimeIdentifier  ("purchaseTime"),
                                                    purchaseTokenIdentifier ("purchaseToken");

                            var orderId          = props[orderIdIdentifier];
                            var appPackageName   = props[packageNameIdentifier];
                            var productId        = props[productIdIdentifier];
                            var purchaseTime     = props[purchaseTimeIdentifier];
                            var purchaseToken    = props[purchaseTokenIdentifier];

                            String purchaseTimeString = Time (purchaseTime.toString().getLargeIntValue()).toString (true, true, true, true);
                            purchases.add ({ orderId.toString(), productId.toString(), appPackageName.toString(), purchaseTimeString, purchaseToken.toString() });
                        }
                    }

                    if (newContinuationToken.get() != 0)
                        getProductsBought (productType, newContinuationToken.get());
                }
            }

            return purchases;
        }

        Pimpl& owner;
        GlobalRef packageName;
        Callback callback;
    };

    //==============================================================================
    class ConsumePurchaseJob : public ThreadPoolJob
    {
    public:
        struct Result
        {
            bool success = false;
            String productIdentifier;
            String statusDescription;
        };

        using Callback = std::function<void(const Result&)>;

        ConsumePurchaseJob (Pimpl& parent,
                            const LocalRef<jstring>& packageNameToUse,
                            const String& productIdentifierToUse,
                            const String& purchaseTokenToUse,
                            const Callback& callbackToUse)
            : ThreadPoolJob ("ConsumePurchaseJob"),
              owner (parent),
              packageName (LocalRef<jobject> (packageNameToUse.get())),
              productIdentifier (productIdentifierToUse),
              purchaseToken (purchaseTokenToUse),
              callback (callbackToUse)
        {}

        ThreadPoolJob::JobStatus runJob() override
        {
            jassert (callback);

            if (owner.checkIsReady())
            {
                auto token = (! purchaseToken.isEmpty() ? purchaseToken : getPurchaseTokenForProductId (productIdentifier, false, 0));

                if (token.isEmpty())
                {
                    if (callback)
                        callback ({ false, productIdentifier, NEEDS_TRANS ("Item not owned") });

                    return jobHasFinished;
                }

                auto responseCode = owner.inAppBillingService.callIntMethod (IInAppBillingService.consumePurchase, 3,
                                                                       (jstring)packageName.get(), javaString (token).get());

                if (callback)
                    callback ({ responseCode == 0, productIdentifier, statusCodeToUserString (responseCode) });
            }
            else
            {
                if (callback)
                    callback ({false, {}, "In-App purchases unavailable"});
            }

            return jobHasFinished;
        }

    private:
        String getPurchaseTokenForProductId (const String productIdToLookFor, bool isSubscription, jstring continuationToken)
        {
            auto productTypeString = javaString (isSubscription ? "subs" : "inapp");
            auto ownedItems = LocalRef<jobject> (owner.inAppBillingService.callObjectMethod (IInAppBillingService.getPurchases, 3,
                                                                                       (jstring) packageName.get(), productTypeString.get(),
                                                                                       continuationToken));

            if (ownedItems.get() != 0)
            {
                auto* env = getEnv();

                auto responseCodeString = javaString ("RESPONSE_CODE");
                auto responseCode = env->CallIntMethod (ownedItems.get(), AndroidBundle.getInt, responseCodeString.get());

                if (responseCode == 0)
                {
                    auto dataListString          = javaString ("INAPP_PURCHASE_DATA_LIST");
                    auto continuationTokenString = javaString ("INAPP_CONTINUATION_TOKEN");

                    auto purchaseDataList     = LocalRef<jobject> (env->CallObjectMethod (ownedItems.get(), AndroidBundle.getStringArrayList, dataListString.get()));
                    auto newContinuationToken = LocalRef<jstring> ((jstring) env->CallObjectMethod (ownedItems.get(), AndroidBundle.getString, continuationTokenString.get()));

                    for (auto i = 0; i < env->CallIntMethod (purchaseDataList.get(), JavaArrayList.size); ++i)
                    {
                        auto purchaseData = juceString ((jstring) (env->CallObjectMethod (purchaseDataList.get(), JavaArrayList.get, i)));

                        var responseData = JSON::parse (purchaseData);

                        if (auto* object = responseData.getDynamicObject())
                        {
                            static const Identifier productIdIdentifier     ("productId"),
                                                    purchaseTokenIdentifier ("purchaseToken");

                            auto& props = object->getProperties();
                            var productId = props[productIdIdentifier];

                            if (productId.toString() == productIdToLookFor)
                                return props[purchaseTokenIdentifier].toString();
                        }
                    }

                    if (newContinuationToken.get() != 0)
                        return getPurchaseTokenForProductId (productIdToLookFor, isSubscription, newContinuationToken.get());
                }
            }

            return {};
        }

        Pimpl& owner;
        GlobalRef packageName;
        const String productIdentifier, purchaseToken;
        Callback callback;
    };

    //==============================================================================
    void handleAsyncUpdate() override
    {
        {
            const ScopedLock lock (getProductsInformationJobResultsLock);

            for (int i = getProductsInformationJobResults.size(); --i >= 0;)
            {
                const auto& result = getProductsInformationJobResults.getReference (i);

                owner.listeners.call ([&] (Listener& l) { l.productsInfoReturned (result); });
                getProductsInformationJobResults.remove (i);
            }
        }

        {
            const ScopedLock lock (getProductsBoughtJobResultsLock);

            for (int i = getProductsBoughtJobResults.size(); --i >= 0;)
            {
                const auto& result = getProductsBoughtJobResults.getReference (i);

                owner.listeners.call ([&] (Listener& l) { l.purchasesListRestored (result.purchases, result.success, result.statusDescription); });
                getProductsBoughtJobResults.remove (i);
            }
        }

        {
            const ScopedLock lock (consumePurchaseJobResultsLock);

            for (int i = consumePurchaseJobResults.size(); --i >= 0;)
            {
                const auto& result = consumePurchaseJobResults.getReference (i);

                owner.listeners.call ([&] (Listener& l) { l.productConsumed (result.productIdentifier, result.success, result.statusDescription); });
                consumePurchaseJobResults.remove (i);
            }
        }
    }

    //==============================================================================
    void inAppPurchaseCompleted (jobject intentData)
    {
        auto* env = getEnv();

        auto inAppPurchaseDataString  = javaString ("INAPP_PURCHASE_DATA");
        auto inAppDataSignatureString = javaString ("INAPP_DATA_SIGNATURE");
        auto responseCodeString       = javaString ("RESPONSE_CODE");

        auto pd  = LocalRef<jstring> ((jstring) env->CallObjectMethod (intentData, AndroidIntent.getStringExtra, inAppPurchaseDataString.get()));
        auto sig = LocalRef<jstring> ((jstring) env->CallObjectMethod (intentData, AndroidIntent.getStringExtra, inAppDataSignatureString.get()));
        auto purchaseDataString  = pd.get()  != 0 ? juceString (pd.get())  : String();
        auto dataSignatureString = sig.get() != 0 ? juceString (sig.get()) : String();

        var responseData = JSON::parse (purchaseDataString);

        auto responseCode = env->CallIntMethod (intentData, AndroidIntent.getIntExtra, responseCodeString.get());
        auto statusCodeUserString = statusCodeToUserString (responseCode);

        if (auto* object = responseData.getDynamicObject())
        {
            auto& props = object->getProperties();

            static const Identifier orderIdIdentifier          ("orderId"),
                                    packageNameIdentifier      ("packageName"),
                                    productIdIdentifier        ("productId"),
                                    purchaseTimeIdentifier     ("purchaseTime"),
                                    purchaseTokenIdentifier    ("purchaseToken"),
                                    developerPayloadIdentifier ("developerPayload");

            var orderId          = props[orderIdIdentifier];
            var packageName      = props[packageNameIdentifier];
            var productId        = props[productIdIdentifier];
            var purchaseTime     = props[purchaseTimeIdentifier];
            var purchaseToken    = props[purchaseTokenIdentifier];
            var developerPayload = props[developerPayloadIdentifier];

            auto purchaseTimeString = Time (purchaseTime.toString().getLargeIntValue())
                                        .toString (true, true, true, true);

            notifyAboutPurchaseResult ({ orderId.toString(), productId.toString(), packageName.toString(),
                                         purchaseTimeString, purchaseToken.toString() },
                                       true, statusCodeUserString);
            return;
        }

        notifyAboutPurchaseResult ({}, false, statusCodeUserString);
    }

    //==============================================================================
    static String statusCodeToUserString (int statusCode)
    {
        switch (statusCode)
        {
            case 0:  return NEEDS_TRANS ("Success");
            case 1:  return NEEDS_TRANS ("Cancelled by user");
            case 2:  return NEEDS_TRANS ("Service unavailable");
            case 3:  return NEEDS_TRANS ("Billing unavailable");
            case 4:  return NEEDS_TRANS ("Item unavailable");
            case 5:  return NEEDS_TRANS ("Internal error");
            case 6:  return NEEDS_TRANS ("Generic error");
            case 7:  return NEEDS_TRANS ("Item already owned");
            case 8:  return NEEDS_TRANS ("Item not owned");
            default: jassertfalse; return NEEDS_TRANS ("Unknown status");
        }
    }

    //==============================================================================
    InAppPurchases& owner;
    GlobalRef inAppBillingService, serviceConnection;
    std::unique_ptr<ThreadPool> threadPool;

    CriticalSection getProductsInformationJobResultsLock,
                    getProductsBoughtJobResultsLock,
                    consumePurchaseJobResultsLock;

    Array<Array<InAppPurchases::Product>> getProductsInformationJobResults;
    Array<GetProductsBoughtJob::Result> getProductsBoughtJobResults;
    Array<ConsumePurchaseJob::Result> consumePurchaseJobResults;
};


//==============================================================================
void juce_inAppPurchaseCompleted (void* intentData)
{
    if (auto* instance = InAppPurchases::getInstance())
        instance->pimpl->inAppPurchaseCompleted (static_cast<jobject> (intentData));
}

} // namespace juce
