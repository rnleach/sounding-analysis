use metfor::Meters;
use optional::Optioned;

/// Station information including location data and identification number.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct StationInfo {
    /// station number, USAF number, eg 727730
    num: Optioned<i32>,
    /// station id, usually a 3 or 4 letter value
    id: Option<String>,
    /// Latitude and longitude.
    location: Option<(f64, f64)>,
    /// Elevation, this may be in model terrain which is not necessarily the same as the real world.
    elevation: Optioned<Meters>,
}

impl StationInfo {
    /// Create a new `StationInfo` object.
    ///
    /// # Arguments
    /// station_num: The USAF station identifier, or None.
    ///
    /// location: The latitude and longitude as a tuple, or None.
    ///
    /// elevation: The elevation of the station **in meters**.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use metfor::{Meters, Feet};
    /// use sounding_analysis::StationInfo;
    /// use optional::{some, none};
    ///
    /// let _stn = StationInfo::new_with_values(12345, None, (45.2,-113.5), Meters(2000.0));
    /// let _stn = StationInfo::new_with_values(12345, None, (45.2,-113.5), Feet(2000.0));
    /// let _stn = StationInfo::new_with_values(12345, None, (45.2,-113.5), some(Meters(2000.0)));
    /// let _stn = StationInfo::new_with_values(12345, None, (45.2,-113.5), some(Feet(2000.0)));
    /// let _stn = StationInfo::new_with_values(12345, None, Some((45.2,-113.5)), Meters(2000.0));
    /// let _stn = StationInfo::new_with_values(12345, None, Some((45.2,-113.5)), Feet(2000.0));
    /// let _stn = StationInfo::new_with_values(12345, None, Some((45.2,-113.5)), some(Meters(2000.0)));
    /// let _stn = StationInfo::new_with_values(12345, None, Some((45.2,-113.5)), some(Feet(2000.0)));
    /// let _stn = StationInfo::new_with_values(Some(12345), None, None, Meters(2000.0));
    /// let _stn = StationInfo::new_with_values(Some(12345), None, None, Feet(2000.0));
    /// let _stn = StationInfo::new_with_values(None, None, (45.2,-113.5), some(Meters(2000.0)));
    /// let _stn = StationInfo::new_with_values(None, None, (45.2,-113.5), some(Feet(2000.0)));
    ///
    /// // Note that lat-lon is an `Option` and not an `Optioned`
    /// let _stn = StationInfo::new_with_values(some(12345), None, None, none::<Feet>());
    /// let _stn = StationInfo::new_with_values(some(12345), None, None, none::<Meters>());
    /// ```
    #[inline]
    pub fn new_with_values<S, T, U, V, W>(
        station_num: T,
        station_id: S,
        location: U,
        elevation: V,
    ) -> Self
    where
        S: Into<Option<String>>,
        T: Into<Optioned<i32>>,
        U: Into<Option<(f64, f64)>>,
        Optioned<W>: From<V>,
        W: optional::Noned + metfor::Length,
        Meters: From<W>,
    {
        let elev: Optioned<W> = Optioned::from(elevation);
        let elev: Optioned<Meters> = elev.map_t(Meters::from);

        StationInfo {
            num: station_num.into(),
            id: station_id.into(),
            location: location.into(),
            elevation: elev,
        }
    }

    /// Create a new object with default values.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use sounding_analysis::StationInfo;
    ///
    /// assert!(StationInfo::new().station_num().is_none());
    /// assert!(StationInfo::new().location().is_none());
    /// assert!(StationInfo::new().elevation().is_none());
    ///
    /// ```
    #[inline]
    pub fn new() -> Self {
        Self::default()
    }

    /// Builder method to add a station number.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use sounding_analysis::StationInfo;
    ///
    /// assert_eq!(StationInfo::new().with_station(12345).station_num().unwrap(), 12345);
    /// assert_eq!(StationInfo::new().with_station(Some(12345)).station_num().unwrap(), 12345);
    ///
    /// ```
    #[inline]
    pub fn with_station<T>(mut self, number: T) -> Self
    where
        Optioned<i32>: From<T>,
    {
        self.num = Optioned::from(number);

        self
    }

    /// Builder method to add a location.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use sounding_analysis::StationInfo;
    ///
    /// assert_eq!(
    ///     StationInfo::new().with_lat_lon((45.0, -116.0)).location().unwrap(), (45.0, -116.0));
    /// assert_eq!(
    ///     StationInfo::new().with_lat_lon(Some((45.0, -116.0)))
    ///         .location()
    ///         .unwrap(),
    ///     (45.0, -116.0));
    ///
    /// ```
    #[inline]
    pub fn with_lat_lon<T>(mut self, coords: T) -> Self
    where
        Option<(f64, f64)>: From<T>,
    {
        self.location = Option::from(coords);
        self
    }

    /// Builder method to add elevation.
    ///
    /// # Examples
    ///```rust
    /// use metfor::{Meters, Feet, Km};
    /// use sounding_analysis::StationInfo;
    /// use optional::{some, none};
    ///
    /// let _info = StationInfo::new().with_elevation(Feet(200.0));
    /// let _info = StationInfo::new().with_elevation(Meters(200.0));
    /// let _info = StationInfo::new().with_elevation(Km(2.0));
    /// let _info = StationInfo::new().with_elevation(some(Feet(200.0)));
    /// let _info = StationInfo::new().with_elevation(some(Meters(200.0)));
    /// let _info = StationInfo::new().with_elevation(some(Km(2.0)));
    /// let _info = StationInfo::new().with_elevation(none::<Feet>());
    /// let _info = StationInfo::new().with_elevation(none::<Meters>());
    /// let _info = StationInfo::new().with_elevation(none::<Km>());
    ///```
    #[inline]
    pub fn with_elevation<T, U>(mut self, elev: T) -> Self
    where
        Optioned<U>: From<T>,
        U: optional::Noned + metfor::Length,
        Meters: From<U>,
    {
        let elevation: Optioned<U> = Optioned::from(elev);
        let elevation: Optioned<Meters> = elevation.map_t(Meters::from);

        self.elevation = elevation;
        self
    }

    /// Builder method to add a station ID. These are usually 3 or 4 alphanumeric codes that may
    /// not be unique to the location like the station number is supposed to be.
    #[inline]
    pub fn with_station_id<T>(mut self, station_id: T) -> Self
    where
        Option<String>: From<T>,
    {
        let id: Option<String> = Option::from(station_id);
        self.id = id;
        self
    }

    /// station number, USAF number, eg 727730
    #[inline]
    pub fn station_num(&self) -> Optioned<i32> {
        self.num
    }

    /// Latitude and longitude.
    #[inline]
    pub fn location(&self) -> Option<(f64, f64)> {
        self.location
    }

    /// Elevation in meters, this may be in model terrain, not necessarily the same as
    /// the real world.
    #[inline]
    pub fn elevation(&self) -> Optioned<Meters> {
        self.elevation
    }

    /// Get the station ID that was used with this station. This is normally a series of 3 or 4
    /// letters. It is not unique to the location like the station number is supposed to be.
    ///
    /// # Examples
    /// ```
    /// use sounding_analysis::StationInfo;
    ///
    /// let info = StationInfo::new().with_station_id("KXLY".to_owned());
    /// assert_eq!(Some("KXLY"), info.station_id());
    /// ```
    #[inline]
    pub fn station_id(&self) -> Option<&str> {
        self.id.as_ref().map(|s| s.as_ref())
    }
}
