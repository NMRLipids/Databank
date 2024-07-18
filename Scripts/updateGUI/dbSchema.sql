-- MySQL dump 10.13  Distrib 8.0.38, for Linux (x86_64)
--
-- Host: localhost    Database: nmrlipid_trajectdb
-- ------------------------------------------------------
-- Server version	8.0.38

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!50503 SET NAMES utf8mb4 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Table structure for table `experiments_FF`
--

DROP TABLE IF EXISTS `experiments_FF`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `experiments_FF` (
  `id` bigint unsigned NOT NULL AUTO_INCREMENT,
  `doi` varchar(255) CHARACTER SET utf8mb3 COLLATE utf8mb3_general_ci NOT NULL,
  `path` varchar(255) CHARACTER SET utf8mb3 COLLATE utf8mb3_general_ci NOT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=39 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `experiments_OP`
--

DROP TABLE IF EXISTS `experiments_OP`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `experiments_OP` (
  `id` bigint unsigned NOT NULL AUTO_INCREMENT,
  `doi` varchar(255) CHARACTER SET utf8mb3 COLLATE utf8mb3_general_ci NOT NULL,
  `path` varchar(255) CHARACTER SET utf8mb3 COLLATE utf8mb3_general_ci NOT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=34 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `forcefields`
--

DROP TABLE IF EXISTS `forcefields`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `forcefields` (
  `id` bigint unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(255) NOT NULL DEFAULT '',
  `date` varchar(255) NOT NULL DEFAULT '',
  `source` varchar(255) DEFAULT NULL,
  `created_at` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `updated_at` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=139 DEFAULT CHARSET=utf8mb3;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `heteromolecules`
--

DROP TABLE IF EXISTS `heteromolecules`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `heteromolecules` (
  `id` bigint unsigned NOT NULL AUTO_INCREMENT,
  `molecule` varchar(255) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci NOT NULL,
  `forcefield_id` bigint unsigned NOT NULL,
  `name` varchar(255) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci NOT NULL,
  `mapping` varchar(255) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci NOT NULL,
  `created_at` timestamp NULL DEFAULT CURRENT_TIMESTAMP,
  `updated_at` timestamp NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`),
  KEY `forcefield_id` (`forcefield_id`),
  CONSTRAINT `heteromolecules_ibfk_1` FOREIGN KEY (`forcefield_id`) REFERENCES `forcefields` (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=37 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ions`
--

DROP TABLE IF EXISTS `ions`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `ions` (
  `id` bigint unsigned NOT NULL AUTO_INCREMENT,
  `molecule` varchar(255) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci NOT NULL,
  `forcefield_id` bigint unsigned NOT NULL,
  `name` varchar(255) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci NOT NULL,
  `mapping` varchar(255) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci NOT NULL,
  `created_at` timestamp NULL DEFAULT NULL,
  `updated_at` timestamp NULL DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `forcefield_id` (`forcefield_id`),
  CONSTRAINT `ions_ibfk_1` FOREIGN KEY (`forcefield_id`) REFERENCES `forcefields` (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=183 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `lipids`
--

DROP TABLE IF EXISTS `lipids`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `lipids` (
  `id` bigint unsigned NOT NULL AUTO_INCREMENT,
  `molecule` varchar(255) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci NOT NULL,
  `forcefield_id` bigint unsigned NOT NULL,
  `name` varchar(255) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci NOT NULL,
  `mapping` varchar(255) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci NOT NULL,
  `created_at` timestamp NULL DEFAULT NULL,
  `updated_at` timestamp NULL DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `forcefield_id` (`forcefield_id`),
  CONSTRAINT `lipids_ibfk_1` FOREIGN KEY (`forcefield_id`) REFERENCES `forcefields` (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=278 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `membranes`
--

DROP TABLE IF EXISTS `membranes`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `membranes` (
  `id` bigint unsigned NOT NULL AUTO_INCREMENT,
  `forcefield_id` bigint unsigned NOT NULL,
  `lipid_names_l1` varchar(255) DEFAULT NULL,
  `lipid_names_l2` varchar(255) DEFAULT NULL,
  `lipid_number_l1` varchar(255) DEFAULT NULL,
  `lipid_number_l2` varchar(255) DEFAULT NULL,
  `geometry` varchar(255) DEFAULT NULL,
  `created_at` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  `updated_at` timestamp NOT NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (`id`),
  KEY `forcefield_id` (`forcefield_id`),
  CONSTRAINT `membranes_ibfk_1` FOREIGN KEY (`forcefield_id`) REFERENCES `forcefields` (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=495 DEFAULT CHARSET=utf8mb3;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `migrations`
--

DROP TABLE IF EXISTS `migrations`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `migrations` (
  `id` int unsigned NOT NULL AUTO_INCREMENT,
  `migration` varchar(255) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci NOT NULL,
  `batch` int NOT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=9 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ranking_global`
--

DROP TABLE IF EXISTS `ranking_global`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `ranking_global` (
  `id` bigint NOT NULL AUTO_INCREMENT,
  `trajectory_id` bigint NOT NULL,
  `ranking_total` bigint NOT NULL,
  `ranking_hg` float NOT NULL,
  `ranking_tails` float NOT NULL,
  `quality_total` float NOT NULL,
  `quality_hg` float NOT NULL,
  `quality_tails` float NOT NULL,
  PRIMARY KEY (`id`) USING BTREE,
  KEY `trajectory_id` (`trajectory_id`)
) ENGINE=MyISAM AUTO_INCREMENT=789 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ranking_heteromolecules`
--

DROP TABLE IF EXISTS `ranking_heteromolecules`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `ranking_heteromolecules` (
  `id` bigint NOT NULL AUTO_INCREMENT,
  `trajectory_id` bigint NOT NULL,
  `molecule_id` bigint NOT NULL,
  `ranking_total` int NOT NULL,
  `quality_total` float NOT NULL,
  PRIMARY KEY (`id`),
  KEY `trajectory_id` (`trajectory_id`),
  KEY `molecule_id` (`molecule_id`)
) ENGINE=InnoDB AUTO_INCREMENT=192 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `ranking_lipids`
--

DROP TABLE IF EXISTS `ranking_lipids`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `ranking_lipids` (
  `id` bigint NOT NULL AUTO_INCREMENT,
  `trajectory_id` bigint NOT NULL,
  `lipid_id` int NOT NULL,
  `ranking_total` int NOT NULL,
  `ranking_hg` int NOT NULL,
  `ranking_sn-1` int NOT NULL,
  `ranking_sn-2` int NOT NULL,
  `quality_total` float NOT NULL,
  `quality_hg` float NOT NULL,
  `quality_sn-1` float NOT NULL,
  `quality_sn-2` float NOT NULL,
  PRIMARY KEY (`id`) USING BTREE,
  KEY `trajectory_id` (`trajectory_id`),
  KEY `lipid_id` (`lipid_id`)
) ENGINE=InnoDB AUTO_INCREMENT=989 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trajectories`
--

DROP TABLE IF EXISTS `trajectories`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `trajectories` (
  `id` bigint unsigned NOT NULL AUTO_INCREMENT,
  `forcefield_id` bigint unsigned NOT NULL,
  `membrane_id` bigint unsigned NOT NULL,
  `git_path` varchar(255) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci NOT NULL,
  `system` varchar(255) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci NOT NULL,
  `author` varchar(255) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci NOT NULL,
  `date` varchar(255) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci NOT NULL,
  `dir_wrk` varchar(255) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci NOT NULL,
  `doi` varchar(255) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci NOT NULL,
  `number_of_atoms` int NOT NULL,
  `preeq_time` varchar(255) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci NOT NULL,
  `publication` varchar(255) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci NOT NULL,
  `temperature` double NOT NULL,
  `software` varchar(255) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci NOT NULL DEFAULT '',
  `trj_size` double NOT NULL,
  `trj_length` double NOT NULL,
  `timeleftout` double NOT NULL,
  `created_at` timestamp NULL DEFAULT NULL,
  `updated_at` timestamp NULL DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `forcefield_id` (`forcefield_id`),
  KEY `membrane_id` (`membrane_id`),
  CONSTRAINT `trajectories_ibfk_1` FOREIGN KEY (`forcefield_id`) REFERENCES `forcefields` (`id`),
  CONSTRAINT `trajectories_ibfk_2` FOREIGN KEY (`membrane_id`) REFERENCES `membranes` (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=814 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trajectories_analysis`
--

DROP TABLE IF EXISTS `trajectories_analysis`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `trajectories_analysis` (
  `id` bigint NOT NULL AUTO_INCREMENT,
  `trajectory_id` bigint unsigned NOT NULL,
  `bilayer_thickness` float DEFAULT NULL,
  `area_per_lipid` float DEFAULT NULL,
  `area_per_lipid_file` varchar(255) NOT NULL,
  `form_factor_file` varchar(255) NOT NULL,
  `quality_total` double NOT NULL,
  `quality_headgroups` double NOT NULL,
  `quality_tails` double NOT NULL,
  `form_factor_quality` double NOT NULL,
  `form_factor_scaling` double NOT NULL,
  `form_factor_experiment` varchar(255) NOT NULL,
  PRIMARY KEY (`id`),
  KEY `trajectory_id` (`trajectory_id`),
  CONSTRAINT `trajectories_analysis_ibfk_1` FOREIGN KEY (`trajectory_id`) REFERENCES `trajectories` (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=789 DEFAULT CHARSET=utf8mb3;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trajectories_analysis_heteromolecules`
--

DROP TABLE IF EXISTS `trajectories_analysis_heteromolecules`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `trajectories_analysis_heteromolecules` (
  `id` bigint unsigned NOT NULL AUTO_INCREMENT,
  `trajectory_id` bigint unsigned NOT NULL,
  `molecule_id` bigint unsigned NOT NULL COMMENT 'd',
  `quality_total` double NOT NULL,
  `quality_hg` double NOT NULL,
  `quality_tails` double NOT NULL,
  `order_parameters_file` varchar(255) NOT NULL,
  `order_parameters_quality` varchar(255) NOT NULL,
  `order_parameters_experiment` varchar(255) NOT NULL,
  PRIMARY KEY (`id`),
  KEY `trajectory_id` (`trajectory_id`),
  CONSTRAINT `trajectories_analysis_heteromolecules_ibfk_1` FOREIGN KEY (`trajectory_id`) REFERENCES `trajectories` (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=192 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trajectories_analysis_ions`
--

DROP TABLE IF EXISTS `trajectories_analysis_ions`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `trajectories_analysis_ions` (
  `id` bigint NOT NULL AUTO_INCREMENT,
  `trajectory_id` bigint NOT NULL,
  `ion_id` bigint NOT NULL,
  `density_file` varchar(255) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `trajectory_id` (`trajectory_id`)
) ENGINE=InnoDB AUTO_INCREMENT=724 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trajectories_analysis_lipids`
--

DROP TABLE IF EXISTS `trajectories_analysis_lipids`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `trajectories_analysis_lipids` (
  `id` bigint NOT NULL AUTO_INCREMENT,
  `trajectory_id` bigint unsigned NOT NULL,
  `lipid_id` bigint unsigned NOT NULL,
  `quality_total` varchar(255) NOT NULL,
  `quality_hg` varchar(255) NOT NULL,
  `quality_sn-1` varchar(255) NOT NULL,
  `quality_sn-2` varchar(255) NOT NULL,
  `order_parameters_file` varchar(255) NOT NULL,
  `order_parameters_quality` varchar(255) NOT NULL,
  `order_parameters_experiment` varchar(255) NOT NULL,
  PRIMARY KEY (`id`),
  KEY `trajectory_id` (`trajectory_id`),
  KEY `lipid_id` (`lipid_id`),
  CONSTRAINT `trajectories_analysis_lipids_ibfk_1` FOREIGN KEY (`trajectory_id`) REFERENCES `trajectories` (`id`),
  CONSTRAINT `trajectories_analysis_lipids_ibfk_2` FOREIGN KEY (`lipid_id`) REFERENCES `lipids` (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=989 DEFAULT CHARSET=utf8mb3;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trajectories_analysis_water`
--

DROP TABLE IF EXISTS `trajectories_analysis_water`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `trajectories_analysis_water` (
  `id` int NOT NULL AUTO_INCREMENT,
  `trajectory_id` int NOT NULL,
  `water_id` int NOT NULL,
  `density_file` varchar(255) CHARACTER SET latin1 COLLATE latin1_swedish_ci DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `id` (`id`),
  KEY `id_2` (`id`),
  KEY `trajectory_id` (`trajectory_id`)
) ENGINE=InnoDB AUTO_INCREMENT=788 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trajectories_experiments_FF`
--

DROP TABLE IF EXISTS `trajectories_experiments_FF`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `trajectories_experiments_FF` (
  `id` bigint unsigned NOT NULL AUTO_INCREMENT,
  `trajectory_id` bigint unsigned NOT NULL,
  `experiment_id` bigint unsigned NOT NULL,
  PRIMARY KEY (`id`),
  KEY `trajectory_id` (`trajectory_id`),
  KEY `experiment_id` (`experiment_id`),
  CONSTRAINT `experiment_id` FOREIGN KEY (`experiment_id`) REFERENCES `experiments_FF` (`id`),
  CONSTRAINT `trajectory_id` FOREIGN KEY (`trajectory_id`) REFERENCES `trajectories` (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=188 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trajectories_experiments_OP`
--

DROP TABLE IF EXISTS `trajectories_experiments_OP`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `trajectories_experiments_OP` (
  `id` bigint unsigned NOT NULL AUTO_INCREMENT,
  `trajectory_id` bigint unsigned NOT NULL,
  `lipid_id` bigint unsigned NOT NULL,
  `experiment_id` bigint unsigned NOT NULL,
  PRIMARY KEY (`id`),
  KEY `trajectory_id` (`trajectory_id`),
  KEY `lipid_id` (`lipid_id`),
  KEY `experiment_id` (`experiment_id`),
  CONSTRAINT `trajectories_experiments_OP_ibfk_1` FOREIGN KEY (`trajectory_id`) REFERENCES `trajectories` (`id`),
  CONSTRAINT `trajectories_experiments_OP_ibfk_2` FOREIGN KEY (`lipid_id`) REFERENCES `lipids` (`id`),
  CONSTRAINT `trajectories_experiments_OP_ibfk_3` FOREIGN KEY (`experiment_id`) REFERENCES `experiments_OP` (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=162 DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trajectories_heteromolecules`
--

DROP TABLE IF EXISTS `trajectories_heteromolecules`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `trajectories_heteromolecules` (
  `id` bigint unsigned NOT NULL AUTO_INCREMENT,
  `trajectory_id` bigint unsigned NOT NULL,
  `molecule_id` bigint unsigned NOT NULL,
  `molecule_name` varchar(255) CHARACTER SET utf8mb3 COLLATE utf8mb3_unicode_ci NOT NULL,
  `leaflet_1` int NOT NULL,
  `leaflet_2` int NOT NULL,
  `created_at` timestamp NULL DEFAULT NULL,
  `updated_at` timestamp NULL DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `analysis_trajectory_id_foreign` (`trajectory_id`),
  KEY `Molecule_ID` (`molecule_id`) USING BTREE,
  CONSTRAINT `trajectories_heteromolecules_ibfk_1` FOREIGN KEY (`trajectory_id`) REFERENCES `trajectories` (`id`),
  CONSTRAINT `trajectories_heteromolecules_ibfk_2` FOREIGN KEY (`molecule_id`) REFERENCES `heteromolecules` (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=192 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trajectories_ions`
--

DROP TABLE IF EXISTS `trajectories_ions`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `trajectories_ions` (
  `id` bigint unsigned NOT NULL AUTO_INCREMENT,
  `trajectory_id` bigint unsigned NOT NULL,
  `ion_id` bigint unsigned NOT NULL,
  `ion_name` varchar(255) CHARACTER SET utf8mb3 COLLATE utf8mb3_unicode_ci NOT NULL DEFAULT 'from ions select name',
  `number` int NOT NULL,
  `created_at` timestamp NULL DEFAULT NULL,
  `updated_at` timestamp NULL DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `analysis_trajectory_id_foreign` (`trajectory_id`),
  KEY `Ion_ID` (`ion_id`) USING BTREE,
  CONSTRAINT `trajectories_ions_ibfk_1` FOREIGN KEY (`trajectory_id`) REFERENCES `trajectories` (`id`),
  CONSTRAINT `trajectories_ions_ibfk_2` FOREIGN KEY (`ion_id`) REFERENCES `ions` (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=724 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trajectories_lipids`
--

DROP TABLE IF EXISTS `trajectories_lipids`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `trajectories_lipids` (
  `id` bigint unsigned NOT NULL AUTO_INCREMENT,
  `trajectory_id` bigint unsigned NOT NULL,
  `lipid_id` bigint unsigned NOT NULL,
  `lipid_name` varchar(255) CHARACTER SET utf8mb3 COLLATE utf8mb3_unicode_ci NOT NULL,
  `leaflet_1` int NOT NULL,
  `leaflet_2` int NOT NULL,
  `created_at` timestamp NULL DEFAULT NULL,
  `updated_at` timestamp NULL DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `analysis_trajectory_id_foreign` (`trajectory_id`),
  KEY `Lipid_ID` (`lipid_id`),
  CONSTRAINT `trajectories_lipids_ibfk_1` FOREIGN KEY (`trajectory_id`) REFERENCES `trajectories` (`id`),
  CONSTRAINT `trajectories_lipids_ibfk_2` FOREIGN KEY (`lipid_id`) REFERENCES `lipids` (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=989 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trajectories_membranes`
--

DROP TABLE IF EXISTS `trajectories_membranes`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `trajectories_membranes` (
  `id` bigint unsigned NOT NULL AUTO_INCREMENT,
  `trajectory_id` bigint unsigned NOT NULL,
  `membrane_id` bigint unsigned NOT NULL,
  `bulk` int DEFAULT NULL,
  `name` varchar(1024) DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `trajectory_id` (`trajectory_id`,`membrane_id`),
  KEY `membrane_id` (`membrane_id`),
  CONSTRAINT `trajectories_membranes_ibfk_1` FOREIGN KEY (`trajectory_id`) REFERENCES `trajectories` (`id`),
  CONSTRAINT `trajectories_membranes_ibfk_2` FOREIGN KEY (`membrane_id`) REFERENCES `membranes` (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=789 DEFAULT CHARSET=utf8mb3;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `trajectories_water`
--

DROP TABLE IF EXISTS `trajectories_water`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `trajectories_water` (
  `id` bigint unsigned NOT NULL AUTO_INCREMENT,
  `trajectory_id` bigint unsigned NOT NULL,
  `water_id` bigint unsigned NOT NULL,
  `water_name` varchar(255) CHARACTER SET utf8mb3 COLLATE utf8mb3_unicode_ci NOT NULL,
  `number` int NOT NULL,
  `created_at` timestamp NULL DEFAULT NULL,
  `updated_at` timestamp NULL DEFAULT NULL,
  PRIMARY KEY (`id`),
  KEY `analysis_trajectory_id_foreign` (`trajectory_id`),
  KEY `Water_ID` (`water_id`) USING BTREE,
  CONSTRAINT `trajectories_water_ibfk_1` FOREIGN KEY (`trajectory_id`) REFERENCES `trajectories` (`id`),
  CONSTRAINT `trajectories_water_ibfk_2` FOREIGN KEY (`water_id`) REFERENCES `water_models` (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=788 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `users`
--

DROP TABLE IF EXISTS `users`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `users` (
  `id` bigint unsigned NOT NULL AUTO_INCREMENT,
  `name` varchar(255) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci NOT NULL,
  `email` varchar(255) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci NOT NULL,
  `password` varchar(255) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci NOT NULL,
  `remember_token` varchar(100) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci DEFAULT NULL,
  `created_at` timestamp NULL DEFAULT NULL,
  `updated_at` timestamp NULL DEFAULT NULL,
  PRIMARY KEY (`id`),
  UNIQUE KEY `usuarios_email_unique` (`email`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Table structure for table `water_models`
--

DROP TABLE IF EXISTS `water_models`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `water_models` (
  `id` bigint unsigned NOT NULL AUTO_INCREMENT,
  `short_name` text CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci NOT NULL,
  `mapping` varchar(255) CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci NOT NULL,
  `created_at` timestamp NULL DEFAULT NULL,
  `updated_at` timestamp NULL DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=19 DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
/*!40101 SET character_set_client = @saved_cs_client */;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;

/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;

-- Dump completed on 2024-07-18 17:27:46
