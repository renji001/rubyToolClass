class Trans
  PI = 3.14159265358979324
  XPI = 3.14159265358979324 * 3000.0 / 180.0

  def trans(x, y)
    cg02 = wgs2cg02(x, y)
    bd_encrypt(cg02[:lat], cg02[:lon])
  end

  # GCJ-02 to BD-09
  def bd_encrypt (gcjLat, gcjLon)
    x = gcjLon.to_f
    y = gcjLat.to_f
    z = Math.sqrt(x * x + y * y) + 0.00002 * Math.sin(y * XPI)
    theta = Math.atan2(y, x) + 0.000003 * Math.cos(x * XPI)
    bdLon = z * Math.cos(theta) + 0.0065
    bdLat = z * Math.sin(theta) + 0.006
    { 'lat': bdLat, 'lon': bdLon }
  end

  # bd09 转 gcj02
  def bd_decrypt(bdLat, bdLon)
    x = bdLon.to_f - 0.0065, y = bdLat.to_f - 0.006
    z = Math.sqrt(x * x + y * y) - 0.00002 * Math.sin(y * XPI)
    theta = Math.atan2(y, x) - 0.000003 * Math.cos(x * XPI)
    gcjLon = z * Math.cos(theta)
    gcjLat = z * Math.sin(theta)
    { 'lat': gcjLat, 'lon': gcjLon }
  end

  # gcj02 转 wgs84
  def gcj_decrypt_exact(gcjLat, gcjLon)
    initDelta = 0.01
    threshold = 0.000000001
    dLat = initDelta, dLon = initDelta
    mLat = gcjLat.to_f - dLat, mLon = gcjLon.to_f - dLon
    pLat = gcjLat.to_f + dLat, pLon = gcjLon.to_f + dLon
    wgsLat, wgsLon, i = 0
    while true do
      wgsLat = (mLat + pLat) / 2
      wgsLon = (mLon + pLon) / 2
      tmp = wgs2cg02(wgsLat, wgsLon)
      dLat = tmp[:lat] - gcjLat.to_f
      dLon = tmp[:lon] - gcjLon.to_f
      if ((dLat.abs < threshold) && (dLon.abs < threshold))
        break
      end
      dLat > 0 ? pLat = wgsLat : mLat = wgsLat
      dLon > 0 ? pLon = wgsLon : mLon = wgsLon
      break if (++i > 10000)
    end
    { 'lat': wgsLat, 'lon': wgsLon }
  end

  # WGS84 to GCJ02
  def wgs2cg02(wgsLat, wgsLon)
    d = delta(wgsLat, wgsLon)
    { 'lat': wgsLat.to_f + d[:lat], 'lon': wgsLon.to_f + d[:lon] }
  end

  # wgs84 转 墨卡托
  def wgs2mercator(wgsLat, wgsLon)
    x = wgsLon * 20037508.34 / 180
    y = Math.log(Math.tan((90.+ wgsLat) * PI / 360)) / (PI / 180)
    y = y * 20037508.34 / 180
    { 'lat': y, 'lon': x }
  end

  # 墨卡托投影转 wgs84
  def mercator_decrypt(mercatorLat, mercatorLon)
    x = mercatorLon / 20037508.34 * 180
    y = mercatorLat / 20037508.34 * 180
    y = 180 / PI * (2 * Math.atan(Math.exp(y * PI / 180)) - PI / 2)
    { 'lat': y, 'lon': x };
  end

  def delta(lat, lon)
    lat = lat.to_f
    lon = lon.to_f
    a = 6378245.0; #  a: 卫星椭球坐标投影到平面地图坐标系的投影因子。
    ee = 0.00669342162296594323; #  ee: 椭球的偏心率。
    dLat = transformLat(lon - 105.0, lat - 35.0)
    dLon = transformLon(lon - 105.0, lat - 35.0)
    radLat = lat / 180.0 * PI
    magic = Math.sin(radLat)
    magic = 1 - ee * magic * magic
    sqrtMagic = Math.sqrt(magic)
    dLat = (dLat * 180.0) / ((a * (1 - ee)) / (magic * sqrtMagic) * PI)
    dLon = (dLon * 180.0) / (a / sqrtMagic * Math.cos(radLat) * PI)
    { 'lat': dLat, 'lon': dLon }
  end

  def out_of_china (lat, lon)
    if lon < 72.004 || lon > 137.8347
      return true
    end
    if lat < 0.8293 || lat > 55.8271
      return true
    end
    false
  end

  def transformLat(x, y)
    x = x.to_f
    y = y.to_f
    ret = -100.0 + 2.0 * x + 3.0 * y + 0.2 * y * y + 0.1 * x * y + 0.2 * Math.sqrt(x.abs)
    ret += (20.0 * Math.sin(6.0 * x * PI) + 20.0 * Math.sin(2.0 * x * PI)) * 2.0 / 3.0
    ret += (20.0 * Math.sin(y * PI) + 40.0 * Math.sin(y / 3.0 * PI)) * 2.0 / 3.0
    ret += (160.0 * Math.sin(y / 12.0 * PI) + 320 * Math.sin(y * PI / 30.0)) * 2.0 / 3.0
    ret
  end

  def transformLon(x, y)
    x = x.to_f
    y = y.to_f
    ret = 300.0 + x + 2.0 * y + 0.1 * x * x + 0.1 * x * y + 0.1 * Math.sqrt(x.abs)
    ret += (20.0 * Math.sin(6.0 * x * PI) + 20.0 * Math.sin(2.0 * x * PI)) * 2.0 / 3.0
    ret += (20.0 * Math.sin(x * PI) + 40.0 * Math.sin(x / 3.0 * PI)) * 2.0 / 3.0
    ret += (150.0 * Math.sin(x / 12.0 * PI) + 300.0 * Math.sin(x / 30.0 * PI)) * 2.0 / 3.0
    ret
  end

end